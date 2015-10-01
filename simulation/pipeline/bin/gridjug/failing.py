# coding: utf-8

from time import sleep
from jug import TaskGenerator


@TaskGenerator
def is_prime(n):
    sleep(0.1)

    if n == 6:
        raise RuntimeError

    for j in range(2, n - 1):
        if (n % j) == 0:
            return False
    return True

primes10 = list(map(is_prime, range(2, 11)))
