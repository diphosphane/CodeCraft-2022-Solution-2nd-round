from typing import List, Tuple
from subprocess import getoutput
import numpy as np

LOCAL = getoutput('uname') == 'Darwin'

def read_demand():
    fname = 'data/demand.csv'
    if not LOCAL: fname = '/' + fname
    f = open(fname, 'r'); data = f.readlines(); f.close()
    client_name = data[0].strip().split(',')[2:]
    client_demand = []; time_label = []  # for each time
    stream_id_2_idx_dict = {}; stream_cnt = 0 # stream_id -> idx
    stream_id_list = []
    prev_time = 'the previous time'
    time_record = {}
    for each in data[1:]:
        d = each.strip().split(',')
        curr_time_label = d[0]
        curr_stream_id = d[1]
        if prev_time != curr_time_label:
            prev_time = curr_time_label
            time_label.append(curr_time_label)
            if time_record:
                client_demand.append(time_record)
            time_record = {}  # id -> arr of demand
        if curr_stream_id not in stream_id_2_idx_dict:
            stram_id_idx = stream_cnt
            stream_id_2_idx_dict[curr_stream_id] = stream_cnt
            stream_id_list.append(curr_stream_id)
            stream_cnt += 1
        else:
            stram_id_idx = stream_id_2_idx_dict[curr_stream_id]
        curr_demand = list(map(int, d[2:]))
        time_record[stram_id_idx] = curr_demand
    client_demand.append(time_record)
    return time_label, client_name, client_demand, stream_id_2_idx_dict, stream_id_list

def read_server_bandwidth():
    fname = 'data/site_bandwidth.csv'
    if not LOCAL: fname = '/' + fname
    server_name = []; server_bandwidth = []
    f = open(fname, 'r'); data = f.readlines(); f.close()
    for each in data[1:]:
        sname, bw4server = each.strip().split(',')
        server_name.append(sname)
        server_bandwidth.append(int(bw4server))
    return server_name, server_bandwidth

def read_qos():
    fname = 'data/qos.csv'
    if not LOCAL: fname = '/' + fname
    f = open(fname, 'r'); data = f.readlines(); f.close()
    cname = data[0].strip().split(',')[1:]
    sname = []; qos_array4server = []
    for each in data[1:]:
        qos_line_split = each.strip().split(',')
        sname.append(qos_line_split[0])
        qos_array4server.append(list(map(int, qos_line_split[1:])))
    return cname, sname, qos_array4server

def read_qos_limit():
    fname = 'data/config.ini'
    if not LOCAL: fname = '/' + fname
    f = open(fname, 'r'); data = f.readlines(); f.close()
    qos_lim = int(data[1].strip().split('=')[1])
    base_cost = int(data[2].strip().split('=')[1])
    return qos_lim, base_cost

if __name__ == '__main__':
    read_demand()
    read_qos()
    read_qos_limit()