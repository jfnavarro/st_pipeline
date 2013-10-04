#!/usr/bin/python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

class jc_wrapper(object):
    """
    Simple wrapper to provide a dict-list __setitem__ method to JobConf
    """
    def __init__(self, jc):
        self.jc = jc
        self.cache = {}

    def __getitem__(self, k):
        if self.cache.has_key(k):
            return self.cache[k]
        else:
            return self.jc.get(k)

    def get(self, k):
        return self[k]

    def getInt(self, k):
        return int(self[k])

    def getFloat(self, k):
        return float(self[k])

    def getBoolean(self, k):
        return bool(self[k])

    def __setitem__(self, k, v):
        self.cache[k] = v

    def hasKey(self, k):
        return self.cache.has_key(k) or self.jc.hasKey(k)

# Warn the user that he's using a deprecated property.
# If new_property is provided, the message will suggest to the user to substitute uses
# of the deprecatedProperty with the new one.
def deprecation_warning(log, deprecated_property, new_property):
    log.warning("Your configuration is using the deprecated property %s", deprecated_property)
    if new_property is None:
        log.warning("You should update your configuration to avoid using it.  See the documentation for details.")
    else:
        log.warning("You should update your configuration to replace it with %s. See the documentation for details.", new_property)
        
# Check whether a deprecated property is used, and if so emit a warning.
# @param job_conf A Pydoop JobConf object.
# @param deprecatedProperty The name of the deprecated property.
# @param new_property The name of the property that replaces the deprecated one, if any,
#        to write a suggestion to the user.
def check_deprecated_prop(job_conf, log, deprecated_property, new_property):
    if job_conf.hasKey(deprecated_property):
        deprecation_warning(log, deprecated_property, new_property)
        
def convert_job_conf(jobconf, deprecation_map, logger):
    wrapper = jc_wrapper(jobconf)

    for new, old in deprecation_map.iteritems():
        check_deprecated_prop(wrapper, logger, old, new)
        if wrapper.hasKey(old):
            wrapper[new] = wrapper[old]
            logger.warning("Using value %s for property %s (value taken from its deprecated equivalent property %s).", wrapper[new], new, old)
    return wrapper