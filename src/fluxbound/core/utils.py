import re

VALID_SBML_ID_REGEX = re.compile(r'^[a-zA-Z_][a-zA-Z0-9_]*$')

def valid_sbml_id(x):
    """Test if a string is a valid SBML identifier."""
    return VALID_SBML_ID_REGEX.match(x) is not None


class AttrDict(dict):
    """A dict subclass that allows attribute-style access to its keys."""

    def __getattr__(self, key: str):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"'AttrDict' has no attribute '{key}'") from None

    def __dir__(self):
        return list(self.keys()) + list(super().__dir__())
    

