import logging
import sys

logger = logging.getLogger()
formatter = logging.Formatter("[%(levelname)s] %(asctime)s %(message)s")

logger.handlers = []
console_stream = logging.StreamHandler(sys.stdout)
console_stream.setFormatter(formatter)

logger.addHandler(console_stream)
