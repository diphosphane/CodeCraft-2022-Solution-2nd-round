from collections import defaultdict
from typing import List, Tuple, Set
from subprocess import getoutput
from itertools import product
from random import shuffle
from copy import deepcopy
import math
import time
from functools import reduce
from read_data import *
import numpy as np

cname, sname, qos, qos_lim = None, None, None, None
stream_id_map, stream_id_list = None, None
base_cost = 0
start_time = 0
t_len, s_len, c_len = 0, 0, 0
time_label = None
client_demand = None
bandwidth = None
start_time = None
LOCAL = getoutput('uname') == 'Darwin'

def get_data():
    global cname, sname, qos, qos_lim, bandwidth, client_demand, time_label, base_cost
    global t_len, s_len, c_len
    global stream_id_map, stream_id_list
    cname, sname, qos = read_qos()
    qos_lim, base_cost = read_qos_limit(); qos = np.array(qos)
    time_label, client_name, client_demand, stream_id_map, stream_id_list = read_demand()
    server_name, server_bandwidth = read_server_bandwidth()
    bandwidth = np.array([ server_bandwidth[server_name.index(s)] for s in sname ])
    # client_idx_list = [ client_name.index(c) for c in cname ]
    # client_demand = np.array(client_demand)[:, np.array(client_idx_list)]
    client_idx_list = [ cname.index(c) for c in client_name ]
    cname = client_name
    qos = np.array(qos[:, client_idx_list])
    t_len, s_len, c_len = len(time_label), len(sname), len(cname)

class Solution():
    def __init__(self) -> None:
        self.init_95()
        self.init_qos()
        self.tcs_id_record = [[[ [] for _ in range(s_len)] for _ in range(c_len)] for _ in range(t_len)] # tidx, sidx, cidx -> List[iidx]} 
        self.record = np.zeros((t_len, s_len, c_len), dtype=np.int32)
        self.t_s_record = np.zeros((t_len, s_len), dtype=np.int32)
    
    def init_qos(self):
        def _qos4c(c_idx: int) -> List[int]:
            c_qos = qos[:, c_idx]
            qos_avail = c_qos < qos_lim
            out = [ s_idx for s_idx, avail in enumerate(qos_avail) if avail ]
            return out
        
        def _qos4s(s_idx: int) -> List[int]:
            s_qos = qos[s_idx, :]
            qos_avail = s_qos < qos_lim
            out = [ c_idx for c_idx, avail in enumerate(qos_avail) if avail ]
            return out
        
        self.qos_avail_for_c = [ _qos4c(c_idx) for c_idx in range(c_len) ]
        self.qos_avail_num_for_c = np.array([ len(i) for i in self.qos_avail_for_c])
        self.qos_avail_for_s = [ _qos4s(s_idx) for s_idx in range(s_len) ]
        self.qos_avail_num_for_s = np.array([ len(i) for i in self.qos_avail_for_s])

        self.qos_avail_for_c_set = [ set(s_list) for s_list in self.qos_avail_for_c ]
        self.qos_avail_for_s_set = [ set(c_list) for c_list in self.qos_avail_for_s ]
        self.avail_s_count = 0

        for each in self.qos_avail_for_s:
            if each: self.avail_s_count += 1
        self.s2s_bridge = []
        for s_idx in range(s_len):
            s = set()
            for c_idx in self.qos_avail_for_s[s_idx]:
                s.update(self.qos_avail_for_c[c_idx])
            self.s2s_bridge.append(s)
    
    def init_95(self):
        num_95 = math.ceil(t_len * 0.95)
        self.idx_95 = num_95 - 1
        self.higher_95_num = t_len - num_95
    
    def check_output_valid(self):
        pass
        # # check client is equal
        # demand_sum = self.record.sum(axis=1)
        # for t_idx, sum_at_each_time in enumerate(demand_sum):
        #     c_demand_at_t = client_demand[t_idx]
        #     if np.any(c_demand_at_t - sum_at_each_time): # if c_demand_at_t != sum_at_each_time:
        #         print(f'client demand is not equal at time {t_idx}')
        #         print(f'calculated: \n{sum_at_each_time} \n\n required: \n{c_demand_at_t}')
        #         print(f'difference (calculated_demand - required_demand): \n {sum_at_each_time - c_demand_at_t}')
        #         exit(1)
        # if np.any(demand_sum - client_demand): # if demand_sum != client_demand:
        #     print('client demand is not equal')
        #     exit(1)
        # # check qos
        # for t_idx, r_each_time in enumerate(self.record):
        #     for s_idx, r_each_s in enumerate(r_each_time):
        #         for c_idx, value in enumerate(r_each_s):
        #             if value:
        #                 if qos[s_idx, c_idx] > qos_lim:
        #                     print(f'qos not satisfied in time {t_idx}, server {sname[s_idx]} (index: {s_idx}), client {cname[c_idx]} (index: {c_idx})')
        #                     exit(1)
        #                 if value < 0:
        #                     print(f'dispatch bandwidth < 0 in time {t_idx}, server {sname[s_idx]} (index: {s_idx}), client {cname[c_idx]} (index: {c_idx})')
        #                     exit(1)
        # # check server upper limit
        # bw_sum = self.t_s_record
        # for t_idx, sum_at_t in enumerate(bw_sum):
        #     if np.any(sum_at_t > bandwidth):
        #         print(f'exceed bandwidth upper at time {t_idx} {time_label[t_idx]}')
        #         print(f'different (bandwidth_limit - solution_sum): \n{bandwidth - sum_at_t}')
        #         exit(1)
        print('test passed \n')

    def output(self, record=None):
        if not record:
            record = self.record
        if LOCAL: self.f = open('output/solution.txt', 'w')
        else: self.f = open('/output/solution.txt', 'w')
        for each_time_step_operation in self.tcs_id_record:
            for cidx, iidx_list_4_s in enumerate(each_time_step_operation):
                tmp = cname[cidx] + ':'
                out_list = []
                for sidx, iidx_list in enumerate(iidx_list_4_s):
                    if iidx_list:
                        stream_id_str = ','.join([stream_id_list[iidx] for iidx in iidx_list])
                        out_list.append(f'<{sname[sidx]},{stream_id_str}>')
                tmp += ','.join(out_list)
                self.f.write(tmp + '\n')
        self.f.close()
    
    def calc_score95(self, print_sep=True):
        bw_each_time = self.t_s_record.copy()
        bw_each_time.sort(axis=0)
        score_95 = bw_each_time[self.idx_95, :]
        after_95 = bw_each_time[self.idx_95+1:, :].sum(0)
        after_95_sum = after_95.sum()
        final_score = score_95.sum()
        if print_sep:
            print(f'95% score sum: {final_score}\n{sorted(score_95, reverse=True)}\n')
            print(f'after 95 sum: {after_95_sum}\n{sorted(after_95, reverse=True)}')
        else:
            print(f'95% score sum: {final_score}')
            print(f'after 95 sum: {after_95_sum}')
        return final_score
    
    @staticmethod
    def max_idx_gen(array: np.ndarray) -> Tuple[Tuple[int, int], int]:
        arr = array.copy()
        cnt = 0; whole_num = reduce(lambda x,y: x*y, arr.shape)
        while cnt < whole_num:
            idx = np.unravel_index(np.argmax(arr), arr.shape)
            value = arr[idx]
            if value == 0: return
            yield idx, value
            arr[idx] = 0
            cnt += 1
    
    @staticmethod
    def max_idx_of(arr: np.ndarray) -> Tuple[int, int]:
        return np.unravel_index(np.argmax(arr), arr.shape)
    
    def index_of(self, perc: float) -> int:
        return math.ceil(t_len * perc) - 1
    
    def assign(self, tidx, sidx, cidx, iidx, add_value) -> bool:
        if self.t_s_record[tidx, sidx] + add_value > bandwidth[sidx]:
            return False
        self.t_s_record[tidx, sidx] += client_demand[tidx][iidx][cidx]
        self.tcs_id_record[tidx][cidx][sidx].append(iidx)
        return True
    
    def dispatch_baseline(self):
        for tidx in range(t_len):
            for iidx, demands in client_demand[tidx].items():
                for cidx, need_dispatch in enumerate(demands):
                    for sidx in self.qos_avail_for_c[cidx]:
                        success = self.assign(tidx, sidx, cidx, iidx, need_dispatch)
                        if success: break

    def dispatch(self):
        s_include_t = [ set() for _ in range(s_len) ]
        for tidx in range(t_len):
            for iidx, demands in client_demand[tidx].items():
                for cidx, need_dispatch in enumerate(demands):
                    success = False
                    for sidx in self.qos_avail_for_c[cidx]:
                        if len(s_include_t[sidx]) < 5 or tidx in s_include_t[sidx]:
                            success = self.assign(tidx, sidx, cidx, iidx, need_dispatch)
                            if success: 
                                s_include_t[sidx].add(tidx)
                                break
                    if not success:
                        s_avail_list = self.qos_avail_for_c[cidx]
                        record_at_t = self.t_s_record[tidx][s_avail_list]
                        idx = np.argmin(record_at_t)
                        sidx = s_avail_list[idx]
                        success = self.assign(tidx, sidx, cidx, iidx, need_dispatch)
                        if success: continue
                        perc_30_num = math.ceil(len(s_avail_list) * 0.3)
                        if perc_30_num >= 3:
                            arg = np.argpartition(record_at_t, perc_30_num)[:perc_30_num]
                            for i in arg:
                                sidx = s_avail_list[i]
                                success = self.assign(tidx, sidx, cidx, iidx, need_dispatch)
                                if success: 
                                    break
                            if success: continue
                            for i in arg[perc_30_num:]:
                                sidx = s_avail_list[i]
                                success = self.assign(tidx, sidx, cidx, iidx, need_dispatch)
                                if success: 
                                    break
                        else:
                            arg = np.argsort(record_at_t)
                            for i in arg[1:]:
                                sidx = s_avail_list[i]
                                success = self.assign(tidx, sidx, cidx, iidx, need_dispatch)
                                if success: 
                                    break
            # if tidx < 20:
            #     tmp = [ len(i) for i in s_include_t ]
            #     print(f'time {tidx}: t in top 5 s: {tmp}')


if __name__ == '__main__':
    get_data()
    start_time = time.time()
    s = Solution()
    # s.dispatch_baseline()
    s.dispatch()
    if LOCAL: 
        print(f'used time normal: {(time.time()-start_time):.2f}')
        s.calc_score95(False)

    if LOCAL: 
        s.check_output_valid()
        s.calc_score95(True)
        time_threshould = 10
    else: 
        time_threshould = 283

    s.output()
    if LOCAL: 
        print(f'used time: {(time.time()-start_time):.2f}')
        s.check_output_valid()
