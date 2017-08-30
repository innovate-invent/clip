#!/usr/bin/env python3
import platform
if platform.python_implementation() == 'PyPy':
    import pypysam
else:
    import pysam
from sys import maxsize, stderr
from CigarIterator import CigarIterator, appendOrInc
import multiprocessing, ctypes, io
import parapysam

default_cost = { # Default costs to use when evaluating alignments
    pysam.CMATCH: lambda x: -x,
    pysam.CEQUAL: lambda x: -x,
    pysam.CDIFF: lambda x: -x,
    pysam.CREF_SKIP: lambda x: -x,
    pysam.CINS: lambda x: 6 + 1*(x-1),
    pysam.CDEL: lambda x: 3 + 1*(x-1)
}

def calculateAlignmentCost(record: pysam.AlignedSegment, start: int = 0, end: int = maxsize, costs = default_cost) -> int:
    """
    Calculates the alignment cost of a record.
    :param record: The record to evaluate
    :param start: The start position in cigar space to begin the calculation.
    :param end: The end position in cigar space to end the calculation.
    :param costs: A dictionary mapping op codes to functions that take a single argument of the operation length and returns the cost.
    :return: The total cost of the alignment over the specified range
    """
    cost = 0
    if start < record.reference_start:
        start = record.reference_start
    #if end > record.reference_end:
    #    end = record.reference_end
    if start > end:
        #No work needs to be done
        return cost
    i = CigarIterator(record)
    i.skipToRefPos(start)
    if not i.clipped: cost += costs[i.op](i.opEnd - i.opPos + 1)
    while not i.clipped and i.refPos+i.opLength-1 <= end:
        cost += costs[i.op](i.opLength)
        if not i.nextOp(): break
    if i.valid and not i.clipped:
        cost += costs[i.op](end - i.refPos + 1)
    return cost

def calculateMappingQuality(record: pysam.AlignedSegment) -> int:
    """
    Calculate the mapping quality of the record.
    Currently only returns the current mapping quality as a new one would obfuscate the success of the aligner.
    :param record: The record to calculate the mapping quality
    :return: The mapping quality of the record.
    """
    # Mapping quality isn't updated as it reflects the quality of the original alignment
    return record.mapping_quality

def trimRecord(record: pysam.AlignedSegment, mate: pysam.AlignedSegment, start: int = 0, end: int = maxsize) -> None:
    """
    Soft clip the ends of a record.
    :param record: The record to soft clip
    :param mate: The records mate if applicable or None
    :param start: The start coordinate in reference space of the region of the record to retain.
    :param end: The end coordinate in reference space of the region to retain.
    :return: None
    """
    if start <= record.reference_start and end >= record.reference_end-1:
        # No work needs to be done
        return

    if start < record.reference_start:
        start = record.reference_start
    ops = []

    # Retain the first match after clipping to avoid setting the new start position to a deletion
    nextMatch = None
    i = CigarIterator(record)
    def checkForFirstMatch():
        nonlocal nextMatch, i
        if nextMatch is None and i.op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF):
            nextMatch = i.refPos

    #Jump past hard clipping and retain op
    hard = i.skipClipped(True)
    if hard: ops += [(pysam.CHARD_CLIP, hard)]
    #Skip to reference starting position
    if start <= end and i.skipToRefPos(start):
        #Soft clip everything before start
        appendOrInc(ops, [pysam.CSOFT_CLIP, i.seqPos])
        dist = end - i.refPos #Calculate reference distance remaining to end

        #Copy in all ops before end
        while dist > i.opRemaining or not i.inRef:
            checkForFirstMatch()
            appendOrInc(ops, [i.op, i.opRemaining])
            if i.inRef: dist -= i.opRemaining
            if not i.nextOp(): break

        #If end within op range, copy in remainder
        if i.valid:
            checkForFirstMatch()
            appendOrInc(ops, [i.op, dist])
            if not i.inSeq: dist = 0
            appendOrInc(ops, [pysam.CSOFT_CLIP, i.record.query_length - i.seqPos - dist])
            # Retain hard clip at end
            if len(i.ops) and i.ops[-1][0] == pysam.CHARD_CLIP:
                appendOrInc(ops, i.ops[-1])
    else:
        #Soft clip entire read
        appendOrInc(ops, [pysam.CSOFT_CLIP, record.query_length])
        #Retain hard clip at end
        if len(i.ops) and i.ops[-1][0] == pysam.CHARD_CLIP:
            appendOrInc(ops, i.ops[-1])

    #Update record
    record.cigartuples = ops
    if nextMatch is not None:
        record.set_tag('OS', record.reference_start)
        record.reference_start = nextMatch
    record.mapping_quality = calculateMappingQuality(record)
    # TODO rewrite MD
    if mate: mate.next_reference_start = record.reference_start

def mergeRecord(fromRecord: pysam.AlignedSegment, toRecord: pysam.AlignedSegment, refStart: int = -1, refEnd: int = maxsize, costs = default_cost) -> None:
    """
    Merge the operations of the overlapping region from one record into another. If there is a mismatch of operations at a position between the
    records then the cheapest alignment will be retained.
    :param fromRecord: The record to copy the operations from.
    :param toRecord: The record to copy the operations to.
    :param refStart: The start position of the region to merge in reference space.
                     Will snap to the beginning of the overlapping region if the coordinate is outside of the overlap.
    :param refEnd: The end position of the region to merge in reference space.
                    Will snap to the end of the overlapping region if the coordinate is outside of the overlap.
    :param costs: The cost map to pass to :func:`calculateAlignmentCost`
    :return: None
    """
    overlapStart = fromRecord.reference_start if fromRecord.reference_start > toRecord.reference_start else toRecord.reference_start
    overlapEnd = (fromRecord.reference_end if fromRecord.reference_end < toRecord.reference_end else toRecord.reference_end) - 1

    if refStart < overlapStart:
        refStart = overlapStart
    if refEnd > overlapEnd:
        refEnd = overlapEnd

    if fromRecord.reference_length == 0 or toRecord.reference_length == 0 or refStart >= refEnd:
        #No work needs to be done
        return

    ops = []
    toItr = CigarIterator(toRecord)

    #Copy in unaffected ops in non-overlapping region
    dist = refStart - toItr.refStart
    if not toItr.nextOp(): return
    while not toItr.inRef or dist > toItr.opLength:
        appendOrInc(ops, list(toItr.opRange))
        if toItr.inRef: dist -= toItr.opLength
        if not toItr.nextOp(): return

    appendOrInc(ops, [toItr.op, dist])
    toItr.step(dist)
    seq = toRecord.query_sequence[:toItr.seqPos]
    qual = list(toRecord.query_qualities[:toItr.seqPos])

    #Init fromRecord iterator
    fromItr = CigarIterator(fromRecord)
    if not fromItr.skipToRefPos(refStart): return

    toOptimal = None

    while toItr.refPos <= refEnd and fromItr.refPos <= refEnd:
        if toItr.op == fromItr.op:
            if toItr.inSeq:
                if toItr.baseQual == fromItr.baseQual:
                    r = toItr if not toItr.matchesRef else fromItr
                    seq += r.seqBase #Keep the variant if possible
                    qual += [r.baseQual if r.matchesRef or toItr.seqBase == fromItr.seqBase else 3] # 3 = 50% probability of either base being correct
                elif toItr.baseQual > fromItr.baseQual:
                    seq += toItr.seqBase
                    qual += [toItr.baseQual]
                else:
                    seq += fromItr.seqBase
                    qual += [fromItr.baseQual]
            appendOrInc(ops, [toItr.op, 1])
        else:
            # TODO Merge insertion without moving mate iterator? to try to better align matched regions between mates
            if toOptimal == None:
                # Dont calculate costs if unnecessary
                toOptimal = calculateAlignmentCost(toRecord, refStart, refEnd, costs) < calculateAlignmentCost(fromRecord, refStart, refEnd, costs)
            if toOptimal:
                if toItr.inSeq:
                    seq += toItr.seqBase
                    qual += [toItr.baseQual]
                appendOrInc(ops, [toItr.op, 1])
            else:
                if fromItr.inSeq:
                    seq += fromItr.seqBase
                    qual += [fromItr.baseQual]
                appendOrInc(ops, [fromItr.op, 1])
        if not fromItr.next():
            while toItr.next():  # Copy remainder of toRecord
                if toItr.inSeq:
                    seq += toItr.seqBase
                    qual += [toItr.baseQual]
                appendOrInc(ops, [toItr.op, 1])
            break
        if not toItr.next():
            break
    while toItr.valid:  # Copy remainder of toRecord
        if toItr.inSeq:
            seq += toItr.seqBase
            qual += [toItr.baseQual]
        appendOrInc(ops, [toItr.op, 1])
        toItr.next()
    toRecord.cigartuples = ops
    toRecord.query_sequence = seq
    toRecord.query_qualities = qual
    # TODO update mapping quality
    # TODO Store tag for which strand the base originated from

class WorkerProcess(parapysam.OrderedWorker):
    """
    Worker process that trims and merges records.
    """
    def __init__(self, alternate=False, clipOnly=False, trimTail=False):
        """
        Constructor.
        :param alternate: Set to True to clip trailing mate rather than the first. Use this to try and avoid strand bias.
        :param clipOnly: Set to True to only clip, skipping the merge step.
        :param trimTail: Set to True to trim the tail of the mate record sticking out past the beginning/end of its mate. Use this to remove possible barcode remnants.
        """
        self.alternate = alternate
        self.clipOnly = clipOnly
        self.trimTail = trimTail
        super().__init__()

    def work(self) -> None:
        """
        Called by super class to begin work.
        :return: None
        """
        try:
            firstRecord = self.receiveRecord()
            secondRecord = self.receiveRecord()
            while firstRecord and secondRecord:
                if firstRecord.reference_start < secondRecord.reference_start:
                    leftRecord = firstRecord
                    rightRecord = secondRecord
                else:
                    rightRecord = firstRecord
                    leftRecord = secondRecord

                if not self.clipOnly:
                    if self.alternate:
                        mergeRecord(rightRecord, leftRecord)
                    else:
                        mergeRecord(leftRecord, rightRecord)

                """
                trimTail:
                    If leftRecord.is_reverse: tail is on the left
                    else: tail is on the right
                    If rightRecord.is_reverse: tail is on the left
                    else: tail is on the right
                    if both or neither reversed: skip trimming tail
                    leftRecord:
                        if tail on left: clip left of leftRecord to beginning of rightRecord
                        if tail on right: clip right of leftRecord to end of rightRecord
                    rightRecord:
                        if tail on left: clip left of rightRecord to the beginning of leftRecord
                        if tail on right: clip right of rightRecord to the end of leftRecord
                clip:
                    if alternate: clip left of rightRecord to end of leftRecord
                    else: clip right of leftRecord to beginning of rightRecord
                """
                leftStart = 0
                leftEnd = maxsize
                rightStart = 0
                rightEnd = maxsize

                if self.trimTail and (leftRecord.is_reverse != rightRecord.is_reverse): # XOR
                    if leftRecord.is_reverse:
                        leftStart = rightRecord.reference_start
                    else:
                        leftEnd = rightRecord.reference_end -1

                    if rightRecord.is_reverse:
                        rightStart = leftRecord.reference_start
                    else:
                        rightEnd = leftRecord.reference_end -1

                if self.alternate:
                    rightStart = leftRecord.reference_end
                else:
                    leftEnd = rightRecord.reference_start

                trimRecord(leftRecord, rightRecord, leftStart, leftEnd)
                trimRecord(rightRecord, leftRecord, rightStart, rightEnd)

                self.sendRecord(firstRecord)
                self.sendRecord(secondRecord)
                firstRecord = self.receiveRecord()
                secondRecord = self.receiveRecord()
        except StopIteration:
            pass

class WriterProcess(parapysam.OrderedWorker):
    """
    Writer process that recieves records from worker processes and writes them out.
    """
    def __init__(self, outFH, outFormat, pool, ordered=False):
        """
        Constructor.
        :param outFH: An open file handle to the output file
        :param outFormat: The mode parameter to pass to pysam to determine the output format. Exclude the 'w', it will be prepended.
                          See http://pysam.readthedocs.io/en/latest/api.html#sam-bam-files
        :param pool: A list of the worker processes to read from.
        :param ordered: Set to True to maintain input order, set to 'c' to maintain coordinate order.
        """
        super().__init__()
        self.pool = pool
        self.ordered = ordered
        self.outFH = outFH
        self.outFormat = outFormat
        self.nextIndex = multiprocessing.RawValue(ctypes.c_long, 0)
        self.bufferedCount = multiprocessing.RawValue(ctypes.c_long, 0)

    def work(self) -> None:
        """
        Called by super class to begin work.
        :return: None
        """
        outFile = pysam.AlignmentFile(self.outFH, 'w' + self.outFormat, header=self.header)
        if self.ordered:
            self.writeOrdered(outFile)
        else:
            self.write(outFile)
        outFile.close()

    def stop(self) -> None:
        super().stop(False)
        for p in self.pool:
            p.stop()
        self.join()

    def writeOrdered(self, outFile: pysam.AlignedSegment) -> None:
        """
        Writes output ordered based on :attr:`self.ordered`
        :param outFile: An open pysam.AlignmentFile instance wrapping the output file
        :return: None
        """
        import heapq
        class HeapNode:
            __slots__ = 'key', 'value'
            def __init__(self, key, value):
                self.key, self.value = key, value
            def __lt__(self, other):
                return self.key < other.key
        writeBuffer = []
        self.nextIndex.value = 0  # type: int
        largestPos = 0
        running = True
        while running:
            running = False
            for p in self.pool + [self]: #type: WorkerProcess
                if p.checkEOF():
                    continue
                running = True
                if not p.pollRecord():
                    continue
                result = p.receiveOrderedRecord()
                if self.nextIndex.value == result[0]:
                    self.nextIndex.value += 1

                if self.ordered == 'c':
                    try:
                        startPos = result[1].get_tag('OS')  # type: int
                    except KeyError:
                        startPos = result[1].reference_start
                    if largestPos < startPos:
                        largestPos = startPos
                    heapq.heappush(writeBuffer, HeapNode(result[1].reference_start, result))
                else:
                    heapq.heappush(writeBuffer, HeapNode(result[0], result))
                self.bufferedCount.value += 1
                while len(writeBuffer) and writeBuffer[0].key == self.nextIndex.value and (self.ordered != 'c' or writeBuffer[1].reference_start < largestPos):
                    result = heapq.heappop(writeBuffer)
                    self.bufferedCount.value -= 1
                    if self.ordered != 'c': self.nextIndex.value += 1
                    outFile.write(result.value[1])

        while len(writeBuffer):
            result = heapq.heappop(writeBuffer)
            self.bufferedCount.value -= 1
            if self.ordered != 'c': self.nextIndex.value += 1
            outFile.write(result.value[1])

    def write(self, outFile: pysam.AlignedSegment) -> None:
        """
        Writes records out as they are completed, irrespective of input order.
        :param outFile: An open pysam.AlignmentFile instance wrapping the output file
        :return: None
        """
        running = True
        while running:
            running = False
            for p in self.pool + [self]: #type: WorkerProcess
                if p.checkEOF():
                    continue
                running = True
                if not p.pollRecord():
                    continue
                self.nextIndex.value, record = p.receiveOrderedRecord()
                outFile.write(record)

def _status(mateCount, bufferedCount, nextIndex, logStream: io.IOBase = stderr):
    import time
    logStream.write("\n") # This will be deleted by the next write
    while True:
        status = "\x1b[F\x1b[2K\rMates buffered: {:>10}\tPending output: {:>10}\tWaiting for read number: {:>10}\n".format(mateCount.value, bufferedCount.value, nextIndex.value)
        logStream.write(status)
        time.sleep(0.2)

def clip(inStream: io.IOBase, outStream: io.IOBase, threads: int = 8, maxTLen: int = 1000, outFormat:str = 'bu', ordered=False, alternate=False, clipOnly=False, trimTail=False, verbose=False, logStream: io.IOBase=stderr) -> None:
    """
    Clips all overlapping records from inStream, writing to outstream.
    :param inStream: An open file object pointing to the input file
    :param outStream: An open file object pointing to the output file
    :param threads: The number of worker subprocesses to use. Will reduce to even number if alternate=True.
    :param maxTLen: The maximum expected template length between two mates. If mates exceed this they will not be clipped.
    :param outFormat: The mode parameter to pass to pysam to determine the output format. Exclude the 'w', it will be prepended.
                          See http://pysam.readthedocs.io/en/latest/api.html#sam-bam-files
    :param ordered: Set to True to maintain input order.
    :param alternate: Set to True to clip trailing mate rather than the first. Use this to try and avoid strand bias.
    :param clipOnly: Set to True to only clip, skipping the merge step.
    :param trimTail: Set to True to trim the tail of the mate record sticking out past the beginning/end of its mate. Use this to remove possible barcode remnants.
    :param verbose: Set to True to output status information.
    :param logStream: An open file object to output status information to.
    :return:
    """
    pool = []
    mateBuffer = {}
    inFile = pysam.AlignmentFile(inStream)
    inFileItr = inFile.fetch(until_eof=True)

    # An even number of threads needs to be used to allow round robin loop to alternate consistently at its boundaries
    if alternate:
        threads -= threads % 2
        if threads == 0:
            threads = 2

    # Begin processing
    for i in range(threads):
        worker = WorkerProcess(alternate and bool(i % 2), clipOnly, trimTail)
        worker.start(inFile.header)
        pool.append(worker)

    writer = WriterProcess(outStream, outFormat, pool, 'c' if (alternate or trimTail) and ordered else ordered)
    writer.start(inFile.header)

    def poolLooper():
        while True:
            for p in pool:
                yield p

    poolLoop = poolLooper()
    mateCount = multiprocessing.RawValue(ctypes.c_long, 0)

    statusProc = multiprocessing.Process(target=_status, args=(mateCount, writer.bufferedCount, writer.nextIndex, logStream))
    if verbose:
        statusProc.start()

    i = 0  # Tracks the order the records were read
    for record in inFileItr:
        # Skip clipping code if the records don't overlap
        if not record.is_unmapped \
                and not record.mate_is_unmapped \
                and not record.is_supplementary \
                and record.is_paired \
                and abs(record.template_length) < maxTLen \
                and record.reference_name == record.next_reference_name:
            secondRecord = mateBuffer.get(record.query_name, None)
            if not secondRecord:
                mateBuffer[record.query_name] = (i, record)
                mateCount.value += 1
            else:
                next(poolLoop).sendMatePair(i, record, *secondRecord)
                del mateBuffer[record.query_name]
                mateCount.value -= 1
        else:
            writer.sendOrderedRecord(i, record)
        i += 1
    inFile.close()

    for _, record in mateBuffer.items():
        mateCount.value -= 1
        writer.sendOrderedRecord(*record)
    writer.stop()
    if verbose:
        statusProc.terminate()
        logStream.write("\x1b[F\x1b[2K\rCompleted.\n")
    outStream.close()

if __name__ == '__main__':
    import getopt, os, errno
    from sys import stdout, stdin, argv
    threads = 2
    maxTLen = 1000
    outFormat = ''
    ordered = False
    alternate = False
    trimTail = False
    clipOnly = False
    verbose = False

    #Command line options
    stderr.write("Clip Overlap v1.0\n")
    ops, paths = getopt.gnu_getopt(argv[1:], 't:m:ho:svabc')
    if not len(ops) and not len(paths):
        stderr.write("Use -h for help.\n")
    else:
        for op, val in ops:
            if op == '-h':
                stderr.write("Clip overlapping reads from SAM/BAM/CRAM file\n"
                             "Use: clip.py [-tmos] [input file path | < infile > outfile] [output file path]\n"
                             "If no paths are given stdin and stdout are used.\n"
                             "-t # Threads to use for processing (Default=1)\n"
                             "-m # Maximum template length guaranteeing no read overlap (Default=1000)\n"
                             "-a Alternate strand being clipped to avoid strand bias (RAM intensive)\n"
                             "-b Trim tails of reads that extend past end of mate. Used to trim barcode remnants."
                             "-c Clip only, do not merge clipped region into mate\n"
                             "-o [sbuc] Output format: s=SAM (Default), b=BAM compressed, bu=BAM uncompressed, c=CRAM\n"
                             "-s Maintain input order (High depth regions may fill RAM), if not set will output in arbitrary order (Minimal RAM)\n"
                             "-v Verbose status output\n")
                exit()
            elif op == '-t':
                threads = int(val or 1)
            elif op == '-m':
                maxTLen = int(val or 1000)
            elif op == '-a':
                alternate = True
            elif op == '-b':
                trimTail = True
            elif op == '-c':
                clipOnly = True
            elif op == '-o' and val != 's':
                outFormat += val
            elif op == '-s':
                ordered = True
            elif op == '-v':
                verbose = True

    if len(paths) > 1 and not os.path.exists(os.path.dirname(paths[1])):
        try:
            os.makedirs(os.path.dirname(paths[1]))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    clip(paths[0] if len(paths) else stdin, open(paths[1], 'wb+') if len(paths) > 1 else stdout, threads, maxTLen, outFormat, ordered, alternate, clipOnly, trimTail, verbose)