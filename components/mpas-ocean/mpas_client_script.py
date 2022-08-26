from smartredis import Client
from smartsim import Experiment
from smartsim.database import Orchestrator
from smartsim.log import get_logger, log_to_file
import time
import os
import numpy

#with open('db_debug.log') as f:
#    lines = f.readlines()
#    for line in lines:
#       if 'SSDB' in line:
#           SSDB = line.split('=')[1]
SSDB=os.environ.get("SSDB")
print(f'SSDB in client script {SSDB}')
client = Client(address=SSDB, cluster=False)
recv_array = 20*numpy.ones(1)
err = client.put_tensor("init_recv", recv_array)
key_found = client.poll_key("init_send", 200, 10000)
print('key_found',key_found)
if key_found:
    dummy_array = client.get_tensor("init_send")
    print('init_send_array=',dummy_array[0])
else:
    print('key not found')
