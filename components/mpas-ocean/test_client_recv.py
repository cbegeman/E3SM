from smartredis import Client
from smartsim import Experiment
from smartsim.database import Orchestrator
from smartsim.log import get_logger, log_to_file
import time
import os

client = Client(address="10.128.4.24:6379", cluster=False)
key_found = client.poll_key("send_array", 200, 10000)
print('key_found',key_found)
if key_found:
    dummy_array = client.get_tensor("send_array")
    print('dummy_array=',dummy_array)
else:
    print('key not found')
