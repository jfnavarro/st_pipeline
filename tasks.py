#!/usr/bin/env python

from invoke import run, task
from invoke.util import log

@task
def test():
    """
    Run the tests
    """
    run('python setup.py test', pty=True)

@task
def clean():
    """
    clean - remove build artifacts
    """
    run('rm -rf build/')
    run('rm -rf dist/')
    run('rm -rf stpipeline.egg-info')
    run('rm -rf py-1.4.26-py2.7.egg')
    run('rm -rf pytest-2.6.4-py2.7.egg')
    run('rm -rf test/tmp/test*')
    run('rm -rf test/out/test*')
    run('rm -rf .eggs/')
    run("find . -name '*.pyc' -delete")
    run("find . -name '*.pyo' -delete")
    run("find . -name '*~' -delete")
    run('find . -name __pycache__ -delete')
    log.info('cleaned up')
    
@task(clean)
def publish():
    """
    Publish to the cheeseshop
    """
    run('python setup.py sdist upload', pty=True)
    run('python setup.py bdist_wheel upload', pty=True)
    log.info('published new release')

