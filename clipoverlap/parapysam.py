import multiprocessing, threading, queue, platform
if platform.python_implementation() == 'PyPy':
    import pypysam
else:
    import pysam

class RecordPipe:
    _flushRecordName = "___%!%DummyRecordShouldNotBeInOutput%!%___"
    def __new__(cls, header):
        r, w = multiprocessing.Pipe(False)
        closeEvent = multiprocessing.Event()
        return RecordOutPipe(r, closeEvent), RecordInPipe(w, closeEvent, header)

class RecordInPipe:
    # This is to deal with the fact that pysam locks up when waiting for input but gets EOF
    _flushRecord = pysam.AlignedSegment()
    _flushRecord.query_name = RecordPipe._flushRecordName
    def __init__(self, pipe, closeEvent, header):
        self._pipe = pipe #type: multiprocessing.Connection
        self._closeEvent = closeEvent
        self._header = header
        self._hts = None

    def close(self):
        if self._closeEvent.is_set():
            return
        self.flush()
        if self._hts:
            self._hts.close()
        self._pipe.close()
        self._closeEvent.set()

    @property
    def closed(self):
        return self._closeEvent.is_set()

    def __del__(self):
        self.close()

    def flush(self):
        if not self._hts:
            self._hts = pysam.AlignmentFile(self._pipe, 'wbu', header=self._header)
        for i in range(30):
            self._hts.write(self._flushRecord)

    def write(self, record):
        if not self._hts:
            self._hts = pysam.AlignmentFile(self._pipe, 'wbu', header=self._header)
            #self.flush()
        self._hts.write(record)

class RecordOutPipe:
    def __init__(self, pipe, closeEvent):
        self._pipe = pipe
        self._closeEvent = closeEvent
        self._hts = None
        self._htsItr = None
        self._nextRecord = queue.Queue()
        self._flushing = multiprocessing.Event()
        self._thread = None

    def _initHTS(self):
        if not self._thread:
            self._thread = threading.Thread(target=self._read)
            self._thread.start()

    def close(self):
        if self._hts:
            self._hts.close()
        self._pipe.close()

    @property
    def closed(self):
        return self._closeEvent.is_set()

    @property
    def eof(self):
        return self.closed and self._flushing.is_set() and self._nextRecord.empty()

    def __del__(self):
        self.close()

    def _read(self):
        if not self._htsItr:
            self._hts = pysam.AlignmentFile(self._pipe.fileno(), check_header=False, check_sq=False)
            self._htsItr = self._hts.fetch(until_eof=True)
        try:
            while not self.eof:
                record = next(self._htsItr)
                while record.query_name == RecordPipe._flushRecordName:
                    self._flushing.set()
                    record = next(self._htsItr)
                self._flushing.clear()
                self._nextRecord.put(record)
        except StopIteration:
            pass

    def poll(self):
        self._initHTS()
        return not self._nextRecord.empty()

    def read(self):
        self._initHTS()
        return self._nextRecord.get()

    def __iter__(self):
        return self

    def __next__(self):
        self._initHTS()
        while True:
            try:
                return self._nextRecord.get(True, 1)
            except queue.Empty:
                if self.closed:
                    self.close()
                    raise StopIteration

class WorkerProcess(multiprocessing.Process):
    def start(self, header):
        self.header = header
        self.recordsInR, self.recordsInW = RecordPipe(header)
        self.recordsOutR, self.recordsOutW = RecordPipe(header)
        self._inPID = multiprocessing.current_process().pid
        super().start()

    def run(self):
        self.work()
        self.recordsOutW.close()
        self.recordsInR.close()

    def stop(self, wait=True):
        self.recordsInW.close()
        self.recordsOutR.close()
        if wait: self.join()

    def __del__(self):
        self.stop()

    def sendRecord(self, record):
        if multiprocessing.current_process().pid == self.pid:
            self.recordsOutW.write(record)
        else:
            self.recordsInW.write(record)

    def receiveRecord(self):
        if multiprocessing.current_process().pid == self.pid:
            return next(self.recordsInR)
        else:
            return next(self.recordsOutR)

    def pollRecord(self):
        if multiprocessing.current_process().pid == self.pid:
            return self.recordsInR.poll()
        else:
            return self.recordsOutR.poll()

    def checkEOF(self):
        if multiprocessing.current_process().pid == self.pid:
            return self.recordsInR.eof
        else:
            return self.recordsOutR.eof


class OrderedWorker(WorkerProcess):
    def start(self, header):
        self.orderPipeR, self.orderPipeW = multiprocessing.Pipe(False)
        super().start(header)

    def stop(self, wait=True):
        self.orderPipeW.close()
        super().stop(wait)

    def sendMatePair(self, index1, record1, index2, record2):
        self.sendOrderedRecord(index1, record1)
        self.sendOrderedRecord(index2, record2)

    def sendOrderedRecord(self, index, record):
        self.sendRecord(record)
        self.orderPipeW.send(index)

    def receiveOrderedRecord(self):
        return (self.orderPipeR.recv(), self.receiveRecord())