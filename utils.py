#Here we define util functions
from multiprocessing import Process, Pipe
from itertools import izip


def conjugate(a, b):
    a.conj = b
    if b:
        b.conj = a


def toggleStrand(s):
    if s == '+':
        return '-'
    return '+'


def reverse_complement(seq):
    return seq.translate(
        '*****************************************************************' +
        'TVGHEFCDIJMLKNOPQYSAUBWXRZ[\]^_`tvghefcdijmlknopqysaubwxrz' +
        '*****************************************************************' +
        '+*******************************************************************'
    )[::-1]


def rc(seq):
    return reverse_complement(seq)


def plot(x, y):
    import matplotlib.pyplot as plt
    plt.plot(x, y, '+')
    plt.show()


def spawn(f):
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun


def parmap(f, X):
    pipe = [Pipe() for x in X]
    proc = [Process(target=spawn(f), args=(c, x))
            for x, (p, c) in izip(X, pipe)]
    [p.start() for p in proc]
    [p.join() for p in proc]
    return [p.recv() for (p, c) in pipe]


# lengths should be sorted in increasing order
def Nx(lens, x=50):
    totlen = sum(lens)
    csum = 0
    for i in lens:
        csum += i
        if csum >= totlen * x / 100:
            return i
