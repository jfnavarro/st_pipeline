#!/usr/bin/python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

class EventMonitor(object):
    """
    Interface to count events and measure time spent through counters.
    """
    def __init__(self): pass

    def start(self, event_name):
        """
        Start timing an event.    Stops timing at the corresponding call to
        "stop" with the same event_name
        """
        abstract

    def stop(self, event_name):
        """
        Stop the timer for event_name started with "start"
        """
        abstract

    class TimingBlock(object):
        def __init__(self, monitor, event_name):
            self.__monitor = monitor
            self.__event_name = event_name

        def __enter__(self):
            self.__monitor.start(self.__event_name)
            return self.__monitor

        def __exit__(self, exception_type, exception_val, exception_tb):
            self.__monitor.stop(self.__event_name)
            return False

    def time_block(self, event_name):
        """
        Use in a with statement to time a block of code.  Example:
          with monitor.time_block("my_event"):
              do_something()
        """
        return type(self).TimingBlock(self, event_name)

    def stop_batch(self, event_name, offset, n):
        abstract

    def count(self, event_name, value=1):
        """
        Increment counter "event_name" by "value" (default 1).
        The counter will be created automatically if it
        doesn't already exist.
        """
        abstract

    def has_counter(self, event_name):
        """
        Check whether a counter for "event_name" has already been created.
        """
        abstract

    def add_counter(self, event_name, display_name=None):
        """
        Explicitly create a counter for "event_name", optionally
        setting a display name for reporting.    If no display_name is
        provided, then event_name will be used instead.

        Since counters can be created on-demand by "count", you usually
        won't need to use this method unless you want to set the display_name.
        """
        abstract

    #################################
    # logging methods
    #################################
    def new_status(self, status_msg):
        """
        Report that the job has entered a new state.
        """
        abstract

    def log_debug(self, *args):
        """
        Log a message if log level is at debug or higher.
        """
        abstract

    def log_info(self, *args):
        """
        Log a message if log level is at info or higher.
        """
        abstract

    def log_warning(self, *args):
        """
        Log a message if log level is at warning or higher.
        """
        abstract

    def log_error(self, *args):
        """
        Log a message if log level is at error or higher.
        """
        abstract

    def log_critical(self, *args):
        """
        Log a message if log level is at critical or higher.
        """
        abstract


class QuietMonitor(EventMonitor):
    """
    Implementation of EventMonitor that doesn't do anything.
    """
    def start(self, event_name): pass
    def stop(self, event_name): pass
    def stop_batch(self, event_name, offset, n): pass
    def count(self, event_name, value=1): pass
    def has_counter(self, event_name): pass
    def add_counter(self, event_name, display_name=None): pass
    def new_status(self, status_msg): pass
    def log_debug(self, *args): pass
    def log_info(self, *args): pass
    def log_warning(self, *args): pass
    def log_error(self, *args): pass
    def log_critical(self, *args): pass
