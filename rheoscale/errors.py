

class RheoscaleError(Exception):
    """Base error for RheoScale."""

class ConfigError(RheoscaleError):
    pass

class DataError(RheoscaleError):
    pass

class ComputationError(RheoscaleError):
    pass

class OutputError(RheoscaleError):
    pass
