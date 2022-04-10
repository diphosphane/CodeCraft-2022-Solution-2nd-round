#ifdef __APPLE__
	bool is_local = true;
	#include <iostream> 
	// #include <cstdio> 
	#include <sstream>
	// #include <fstream> 
	#include <algorithm> 
	// #include <cmath> 
	#include <vector>
	// #include <queue> 
	// #include <string> 
	// #include <cstring> 
	#include <map> 
	// #include <stack> 
	#include <set> 
	#include <unordered_set>
#else
	bool is_local = false;
	#include<bits/stdc++.h>
#endif
#include <math.h>
#include <ctime>
#include <float.h>

using namespace std;

///////////////////////////////// global variable definition ///////////////////////////////////////

// input variable
int t_len, s_len, c_len;
vector<string> time_label, cname, sname;
vector<vector<int>> qos;
vector<int> bandwidth;
vector<map<int, vector<int>>> client_demand;
int qos_lim, base_cost;

// intermediate variable
vector<vector<int>> qos4s, qos4c;
vector<int> qos4s_num, qos4c_num;
vector<string> stream_id_vec;
map<string, int> sid_map;
int higher_95_num = 0;
int num_90, num_95, idx_90, idx_95 = 0;
vector<unordered_set<int>> s_include_t;  // temp var
vector<bool> s_used;
vector<int> demand_each_t;
vector<vector<int>> tc_demand;
vector<int> c_demand_one_line;
vector<vector<int>> t_include_sid;
vector<int> barrier;
vector<int> higher_check_num;

// output variable
int ts_record[8928][135];
vector<int> tcs_id_record[8928][35][135];

/////////////////////////////////////////// util //////////////////////////////////////////////////////

inline vector<string> str_split(string &s, char sep) {
	vector<string> out;
	s.erase(s.find_last_not_of("\r")+1);
	stringstream ss(s);
	string item;
	while (getline(ss, item, sep)) {
		out.push_back(item);
	}
	return out;
}

inline int str2num(string &s){
	// stringstream ss(s);
	// int num;
	// ss >> num;
	// return num;
	return stoi(s);
}

template<typename T>
int index_of(vector<T> &vec, T item){
	auto ret = std::find(vec.begin(), vec.end(), item);
	if (ret != vec.end()) return ret - vec.begin();
	return -1;
}

template<class ForwardIter>
inline size_t argmin(ForwardIter first, ForwardIter last){
	return distance(first, min_element(first, last));
}

template<class ForwardIter>
inline size_t argmax(ForwardIter first, ForwardIter last){
	return distance(first, max_element(first, last));
}

template<typename T>
int argmax2(vector<T> &v1, vector<T> &v2){
	T max_value = v1[argmax(v1.begin(), v1.end())];
	vector<int> idx(v1.size());
	for (int i = 0; i < idx.size(); i++) { idx[i] = i; }
	sort(idx.begin(), idx.end(), [&v1, &v2, max_value](int i, int j) -> bool {
		if (v1[i] == max_value && v1[j] < max_value) {
			return true;
		} else if (v1[j] == max_value && v1[i] < max_value) {
			return false;
		}
		return v2[i] > v2[j];
	});
	return idx[0];
}

template<typename T>
vector<int> argsort(vector<T> &v, bool reverse) {
	vector<int> idx(v.size());
	for (int i = 0; i < idx.size(); i++) { idx[i] = i; }
	if (reverse) {
		sort(idx.begin(), idx.end(), [&v](int i, int j) -> bool {
			return v[i] > v[j];
		});
	} else {
		sort(idx.begin(), idx.end(), [&v](int i, int j) -> bool {
			return v[i] < v[j];
		});
	}
	return idx;
}

template<typename T>
vector<int> argsort2(vector<T> &v1, vector<T> &v2, bool reverse) {
	vector<int> idx(v1.size());
	for (int i = 0; i < idx.size(); i++) { idx[i] = i; }
	if (reverse) {
		sort(idx.begin(), idx.end(), [&v1, &v2](int i, int j) -> bool {
			if (v1[i] != v1[j]) {
				return v1[i] > v1[j];
			} else {
				return v2[i] > v2[j];
			}
			
		});
	} else {
		sort(idx.begin(), idx.end(), [&v1, &v2](int i, int j) -> bool {
			if (v1[i] != v1[j]) {
				return v1[i] < v1[j];
			} else {
				return v2[i] < v2[j];
			}
			// return v[i] < v[j];
		});
	}
	return idx;
}

template<typename T>
vector<int> argpartialsort(vector<T> &v, int top_n, bool reverse) {
	vector<int> idx(v.size());
	for (int i = 0; i < idx.size(); i++) { idx[i] = i; }
	if (reverse) {
		nth_element(idx.begin(), idx.begin()+top_n, idx.end(), [&v](int i, int j) -> bool {
				return v[i] > v[j];
		});
		partial_sort(idx.begin(), idx.begin()+top_n, idx.end(), [&v](int i, int j) -> bool {
				return v[i] > v[j];
		});
	} else {
		nth_element(idx.begin(), idx.begin()+top_n, idx.end(), [&v](int i, int j) -> bool {
				return v[i] < v[j];
		});
		partial_sort(idx.begin(), idx.begin()+top_n, idx.end(), [&v](int i, int j) -> bool {
				return v[i] < v[j];
		});
	}
	return idx;
}

template<typename T>
vector<int> argpartition(vector<T> &v, int top_n, bool reverse) {
	vector<int> idx(v.size());
	for (int i = 0; i < idx.size(); i++) { idx[i] = i; }
	if (reverse) {
		nth_element(idx.begin(), idx.begin()+top_n, idx.end(), [&v](int i, int j) -> bool {
				return v[i] > v[j];
		});
	} else {
		nth_element(idx.begin(), idx.begin()+top_n, idx.end(), [&v](int i, int j) -> bool {
				return v[i] < v[j];
		});
	}
	return idx;
}

template<typename T>
int nth_idx(vector<T> &v, int n, bool reverse) {
	vector<int> idx(v.size());
	for (int i = 0; i < idx.size(); i++) { idx[i] = i; }
	if (reverse) {
		nth_element(idx.begin(), idx.begin()+n, idx.end(), [&v](int i, int j) -> bool {
				return v[i] > v[j];
		});
	} else {
		nth_element(idx.begin(), idx.begin()+n, idx.end(), [&v](int i, int j) -> bool {
				return v[i] < v[j];
		});
	}
	return idx[n];
}

//////////////////////////////////////////////// basic part ///////////////////////////////////////////////

void read_qos(){
	if (is_local){
		freopen("data/qos.csv", "r", stdin);
	} else {
		freopen("/data/qos.csv", "r", stdin);
	}
	string s;
	getline(cin, s);
	vector<string> sep_str_vec = str_split(s, ',');
	cname.assign(sep_str_vec.size() - 1, "");
	for (int i=1; i<sep_str_vec.size(); i++) {
		cname[i-1] = sep_str_vec[i];
	}
	c_len = cname.size();
	while (getline(cin, s)) {
		vector<string> sep_str_vec = str_split(s, ',');
		vector<int> qos_new(c_len, 0);
		sname.push_back(sep_str_vec[0]);
		for (int i=1; i<sep_str_vec.size(); i++) {
			qos_new[i-1] = str2num(sep_str_vec[i]);
		}
		qos.push_back(qos_new);
	}
	fclose(stdin);
	cin.clear();
}

void read_config(){
	if (is_local){
		freopen("data/config.ini", "r", stdin);
	} else {
		freopen("/data/config.ini", "r", stdin);
	}
	string s;
	getline(cin, s);
	getline(cin, s);
	qos_lim = str2num(str_split(s, '=')[1]);
	getline(cin, s);
	base_cost = str2num(str_split(s, '=')[1]);
	fclose(stdin);
	cin.clear();
}

void read_bandwidth(){
	if (is_local){
		freopen("data/site_bandwidth.csv", "r", stdin);
	} else {
		freopen("/data/site_bandwidth.csv", "r", stdin);
	}
	vector<string> splitted;
	s_len = sname.size();
	bandwidth.assign(s_len, 0);
	string s;
	int index;
	getline(cin, s);
	while (getline(cin, s)){
		splitted = str_split(s, ',');
		index = index_of(sname, splitted[0]);
		bandwidth[index] = str2num(splitted[1]);
		// bandwidth[index] = stoi(splitted[1]);
	}
	fclose(stdin);
	cin.clear();
}

void re_order(vector<string> &new_cname){
	int index;
	vector<int> vec(c_len, 0), idx_vec(c_len, 0);
	for (int i=0; i < c_len; i++){
		index = index_of(cname, new_cname[i]);
		idx_vec[i] = index;
	}
	for (int i=0; i < s_len; i++){
		vec.assign(c_len, 0);
		for (int j=0; j < c_len; j++){
			vec[j] = qos[i][idx_vec[j]];
		}
		for (int j=0; j < c_len; j++){
			qos[i][j] = vec[j];
		}
	}
	cname = new_cname;
}

void init_data(){
	t_len = time_label.size();
	s_len = sname.size();
	c_len = cname.size();
	s_used.assign(s_len, false);
	for (int i=0; i < s_len; i++){
		unordered_set<int> t_set;
		s_include_t.push_back(t_set);
	}
	num_95 = ceil(t_len * 0.95);
	num_90 = ceil(t_len * 0.9);
	idx_95 = num_95 - 1;
	idx_90 = num_90 - 1;
	higher_95_num = t_len - num_95; 
	barrier.assign(s_len, base_cost);
	higher_check_num.assign(s_len, t_len - num_95);
}

void read_demand(){
	if (is_local){
		freopen("data/demand.csv", "r", stdin);
	} else {
		freopen("/data/demand.csv", "r", stdin);
	}
	string s;
	getline(cin, s);
	vector<string> splitted = str_split(s, ',');
	vector<string> client_name(c_len, "  ");
	for (int i=2; i<splitted.size(); i++){
		client_name[i-2] = splitted[i];
	}
	re_order(client_name);
	cname = client_name;
	string prev_time(" "), time, stream_id;
	vector<int> curr_demand(c_len, 0), c_demand_arr(c_len, 0), sid_idx_arr;
	map<int, vector<int>> timeRecord4id;
	bool first_line = true;
	int num, sid_idx, sid_cnt = 0, demand_accu = 0;
	while (getline(cin, s)){
		splitted = str_split(s, ',');
		time = splitted[0]; stream_id = splitted[1];
		// get and set the sid_idx
		auto find_result = sid_map.find(stream_id);
		if (find_result == sid_map.end()) { // not found
			sid_idx = sid_cnt;
			sid_map[stream_id] = sid_cnt;
			stream_id_vec.push_back(stream_id);
			sid_cnt += 1;
		} else {  // found
			sid_idx = find_result->second;
		}
		// deal with new time
		if (time != prev_time) {
			if (first_line) {
				first_line = false;
			} else {
				client_demand.push_back(timeRecord4id);
				timeRecord4id.clear();
				demand_each_t.push_back(demand_accu);
				demand_accu = 0;
				t_include_sid.push_back(sid_idx_arr);
				sid_idx_arr.clear();
				tc_demand.push_back(c_demand_arr);
				c_demand_one_line.insert(c_demand_one_line.end(), c_demand_arr.begin(), c_demand_arr.end());
				c_demand_arr.assign(c_len, 0);
			}
			time_label.push_back(time);
			prev_time = time;
		}
		// demand convert
		curr_demand.assign(c_len, 0);
		for (int i=2; i<splitted.size(); i++){
			// curr_demand[i-2] = str2num(splitted[i]);
			num = str2num(splitted[i]);
			curr_demand[i-2] = num;
			demand_accu += num;
			c_demand_arr[i-2] += num;
		}
		// store data
		timeRecord4id[sid_idx] = curr_demand;
		sid_idx_arr.push_back(sid_idx);
	}
	client_demand.push_back(timeRecord4id);
	demand_each_t.push_back(demand_accu);
	t_include_sid.push_back(sid_idx_arr);
	tc_demand.push_back(c_demand_arr);
	c_demand_one_line.insert(c_demand_one_line.end(), c_demand_arr.begin(), c_demand_arr.end());
	fclose(stdin);
	cin.clear();
}

void get_data(){
	read_config();
	read_qos();
	read_bandwidth();
	read_demand();
	init_data();
}

void analyse_qos(){
	vector<int> vec;
	for (int i=0; i < c_len; i++){
		vec.clear();
		for (int j=0; j < s_len; j ++) {
			if (qos[j][i] < qos_lim) vec.push_back(j);
		}
		qos4c.push_back(vec);
		qos4c_num.push_back(vec.size());
	}
	for (int i=0; i < s_len; i++){
		vec.clear();
		for (int j=0; j < c_len; j ++) {
			if (qos[i][j] < qos_lim) vec.push_back(j);
		}
		qos4s.push_back(vec);
		qos4s_num.push_back(vec.size());
	}
}

float calc_cost(int s, int w){
	float tmp = pow((float)(w - base_cost), 2) / bandwidth[s] + w;
	return tmp;
}

void output_result(){
	if (is_local){
		freopen("output/solution.txt", "w", stdout);
	} else {
		freopen("/output/solution.txt", "w", stdout);
	}
	// select 90%
	vector<int> value_at_90(s_len, 0), value_at_95(s_len, 0), rec4s(t_len, 0);
	vector<vector<int>> st_record;
	int v90, v95;
	float cost_diff;
	vector<int> cost_diff_arr;
	int num_90 = ceil(t_len * 0.9);
	int num_95 = ceil(t_len * 0.95);
	for (int s = 0; s < s_len; s++) {
		rec4s.assign(t_len, 0);
		for (int t = 0; t < t_len; t++) {
			rec4s[t] = ts_record[t][s];
		}
		st_record.push_back(rec4s);
		v90 = ts_record[nth_idx(rec4s, num_90, false)][s];  // may not need to -1
		v95 = ts_record[nth_idx(rec4s, num_95, false)][s];  // may not need to -1
		cost_diff = calc_cost(s, v95) - calc_cost(s, v90);
		cost_diff_arr.push_back(cost_diff);
		// value_at_90[s] = ts_record[nth_idx(rec4s, num_90, false)-1][s];  // may not need to -1
		// value_at_95[s] = ts_record[nth_idx(rec4s, num_95, false)-1][s];  // may not need to -1
	}
	int can_90_s_num = 0;
	if (is_local) {
		can_90_s_num = 5;
	} else {
		can_90_s_num = 10;
	}
	vector<int> s_idx_arr = argpartition(cost_diff_arr, can_90_s_num, true);
	string tmp;
	int s_idx;
	for (int i = 0; i < can_90_s_num; i++) {
		s_idx = s_idx_arr[i];
		tmp += sname[s_idx] + ',';
	}
	tmp = tmp.substr(0, tmp.length() - 1);
	cout << tmp << endl;
	// output result
	vector<int> sid_idx_arr;
	tmp = "";
	for (int t=0; t < t_len; t++){
		for (int c=0; c < c_len; c++){
			cout << cname[c] << ":";
			tmp = "";
			for (int s=0; s < s_len; s++){
				sid_idx_arr = tcs_id_record[t][c][s];
				if (sid_idx_arr.size() == 0) continue;
				tmp += "<" + sname[s];
				for (auto sid_idx: sid_idx_arr) tmp += "," + stream_id_vec[sid_idx]; //ss4s << "," << stream_id_vec[sid_idx];
				tmp += ">,";
			}
			tmp = tmp.substr(0, tmp.length() - 1);
			cout << tmp << endl;
		}
	}
	cout.clear();
	fclose(stdout);
	freopen("/dev/tty","w",stdout);
}

void write_demand(){
	int cnt = 0;
	for (auto &i: cname) {
		cout << i << "\t";
	}
	cout << endl;
	for (auto &element: client_demand){
		cout << time_label[cnt] << endl;
		for (auto &pair: element) {
			cout << "stream_id ( " << pair.first << " ): ";
			for (auto i: pair.second) cout << i << " ";
			cout << endl;
		}
		if (cnt >= 3) break;
		cout << endl << endl;
		cnt++;
	}
}

///////////////////////////////////////////////// solution 1 baseline ///////////////////////////////////////////////////

bool assign(int tidx, int sidx, int cidx, int sid_idx, int add_value){
	if (ts_record[tidx][sidx] + add_value > bandwidth[sidx]) return false;
	ts_record[tidx][sidx] += add_value;
	tcs_id_record[tidx][cidx][sidx].push_back(sid_idx);
	return true;
}

void dispatch_baseline(){
	int sid_idx;
	vector<int> demand4id;
	int need_dispatch;
	bool success;
	for (int tidx=0; tidx < t_len; tidx++){
		for (auto &demand: client_demand[tidx]){
			sid_idx = demand.first;
			demand4id = demand.second;
			for (int cidx=0; cidx < c_len; cidx++){
				need_dispatch = demand4id[cidx];
				if (need_dispatch == 0) continue;
				for (auto sidx: qos4c[cidx]) {
					success = assign(tidx, sidx, cidx, sid_idx, need_dispatch);
					if (success) break;
				}
			}
		}
	}
}

//////////////////////////////////////////// solution 2 dymanic open server /////////////////////////////////////////////

void free_dispatch(int sid_idx, vector<int> *demand4id, int t, int s) {
	int need_dispatch;
	bool success;
	unordered_set<int> *t_set;
	for (auto c: qos4s[s]){  // free to dispatch  TODO: can be a single function
		need_dispatch = (*demand4id)[c];
		if (need_dispatch == 0) continue;
		t_set = &s_include_t[s];
		// dispatch to include_t = 1 ~ 5
		if ((t_set->size() < higher_95_num) || t_set->find(t) != t_set->end()){
			success = assign(t, s, c, sid_idx, need_dispatch);
			if (success) {
				t_set->insert(t);
				(*demand4id)[c] = 0;  // TODO: may need to delete
			}
		} else { // dispatch to < base_cost  if include_t != 0
			if (need_dispatch + ts_record[t][s] <= barrier[s]) {
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					(*demand4id)[c] = 0;  // TODO: may need to delete
				} else cout << "Fail at ddd";
			}
		}
	}
}

void dispatch_one_stream_id(int sid_idx, vector<int> *demand4id, int t){
	int need_dispatch, min_idx;
	bool success;
	unordered_set<int> *t_set;
	//////////// free to dispatch  TODO: can be a single function /////////////
	for (int c=0; c < c_len; c++){   
		need_dispatch = (*demand4id)[c];
		if (need_dispatch == 0) continue;
		for (auto s: qos4c[c]) {
			if (! s_used[s]) continue;  // filter not used s
			t_set = &s_include_t[s];
			// dispatch to include_t = 1 ~ 5
			if ((t_set->size() < higher_95_num) || t_set->find(t) != t_set->end()){  // can dispatch to 5%
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					t_set->insert(t);
					(*demand4id)[c] = 0;  // TODO: may need to delete
					break;
				}
			} else { // dispatch to < barrier(init is base_cost)  if s is used
				if (need_dispatch + ts_record[t][s] <= barrier[s]) {
					success = assign(t, s, c, sid_idx, need_dispatch);
					if (success) {
						(*demand4id)[c] = 0;  // TODO: may need to delete
						break;
					}
				}
			}
		}
	}
	//////////////////////// pay to dispatch /////////////////////////////////
	int over_flow, cost;
	for (int c=0; c < c_len; c++){
		need_dispatch = (*demand4id)[c];
		if (need_dispatch == 0) continue;
		// cost of dispatch over base cost 
		vector<float> cost_each_s;
		for (auto s: qos4c[c]){  // calculate over flow cost
			if (s_used[s]) {
				t_set = &s_include_t[s];
				if (t_set->find(t) == t_set->end()){  // in 95%, not in 5%
					int sum_here = ts_record[t][s] + need_dispatch;
					if (sum_here > bandwidth[s]) {
						cost_each_s.push_back(FLT_MAX);  // don't use it
					} else {  
						if (ts_record[t][s] < base_cost) {  // only calc the cost that over flow the base cost
							over_flow = sum_here - base_cost;
							cost_each_s.push_back(pow(over_flow, 2) / bandwidth[s] + (sum_here - base_cost));
						} else {  // already >= base_cost
							// cost = dW * ( 1 + (W'+W-2V)/c)  // W: cost at 95%, W': new cost, V: base cost, c: bandwidth
							cost = need_dispatch * (1 + (sum_here + ts_record[t][s] - 2 * base_cost) / bandwidth[s]);
							cost_each_s.push_back(cost);
						}  // WARNING: add up may < barrier, but previous step will prevent this situation
					} 
				} else {  // ???? if in 5%, will be dispatch in previous for_loop
					cost_each_s.push_back(FLT_MAX);
					success = assign(t, s, c, sid_idx, need_dispatch);
					if (success) {
						(*demand4id)[c] = 0;  // TODO: may need to delete
						cout << "WARNING: dispatch at 5%";
					}
					// TODO: try to add value
				}  
			} else {  // not used server, can as a candidate
				cost_each_s.push_back(FLT_MAX / 2);
			}
		}
		min_idx = argmin(cost_each_s.begin(), cost_each_s.end());
		///////////////// over flow current server
		if (cost_each_s[min_idx] < base_cost){  
			int s = qos4c[c][min_idx];
			success = assign(t, s, c, sid_idx, need_dispatch);
			if (success) {
				(*demand4id)[c] = 0;  // TODO: may have BUG
				barrier[s] = max(barrier[s], ts_record[t][s]);
				// TODO: move other res to this server
			} else cout << "Fail at aaa";
		} else {  
		/////////////////// open new server  // TODO: check whether the bandwidth is enough
			vector<int> cand_bw, avail_num;  // choose a new s: 1. bandwidth 2. connectivity
			for (auto s: qos4c[c]){
				if (s_used[s]) cand_bw.push_back(0);
				else {
					cand_bw.push_back(bandwidth[s]);
				}
				avail_num.push_back(qos4s_num[s]);
			}
			// int idx = argmax(cand_bw.begin(), cand_bw.end());
			int idx = argmax2(cand_bw, avail_num);
			int s = qos4c[c][idx];
			if (s_used[s]) {  // s is used 
				s = qos4c[c][min_idx];  // need to overflow the server
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					(*demand4id)[c] = 0;
					barrier[s] = max(barrier[s], ts_record[t][s]);
				// TODO: move other res to this server
				} cout << "Fail at bbb";
			} else {  // open a new s
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					(*demand4id)[c] = 0;
					s_used[s] = true;
					t_set = &s_include_t[s];
					t_set->insert(t);
				} else cout << "Fail at ccc ";  // bandwidth of new s is not enough
				free_dispatch(sid_idx, demand4id, t, s);
			}
		}
	}
}

void dispatch_top5_opt(){
	int sid_idx;
	vector<int> *demand4id;
	vector<int> idxs = argsort(demand_each_t, true);
	for (auto &tidx : idxs) {
	// for (int tidx=0; tidx < t_len; tidx++){
		for (auto &demand: client_demand[tidx]){
			sid_idx = demand.first;
			demand4id = &demand.second;
			dispatch_one_stream_id(sid_idx, demand4id, tidx);
		}
	}
}

///////////////////////////////////////////// solution 3 tc matrix //////////////////////////////////////////////////////

void free_dispatch_a_time(int t, int s, int c) {
	int need_dispatch;
	bool success;
	unordered_set<int> *t_set;
	for (auto &sid_idx : t_include_sid[t]) {
		// need_dispatch = (*demand4id)[c];
		need_dispatch = client_demand[t][sid_idx][c];
		if (need_dispatch == 0) continue;
		t_set = &s_include_t[s];
		// dispatch to include_t = 1 ~ 5
		if ((t_set->size() < higher_95_num) || t_set->find(t) != t_set->end()){
			success = assign(t, s, c, sid_idx, need_dispatch);
			if (success) {
				t_set->insert(t);
				// (*demand4id)[c] = 0;  // TODO: may need to delete
				client_demand[t][sid_idx][c] = 0;
			}
		} else { // dispatch to < base_cost  if include_t != 0
			if (need_dispatch + ts_record[t][s] <= barrier[s]) {
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					// (*demand4id)[c] = 0;  // TODO: may need to delete
					client_demand[t][sid_idx][c] = 0;
				} else cout << "Fail at ddd";
			}
		}
	}
	
}

void free_dispatch_t_c(int t, int c) {
	int need_dispatch, min_idx; //, sid_idx;
	bool success;
	vector<int> *demand4id;
	unordered_set<int> *t_set;
	//////////// free to dispatch  TODO: can be a single function /////////////
	for (auto &sid_idx : t_include_sid[t]) {
		demand4id = &client_demand[t][sid_idx];
		need_dispatch = (*demand4id)[c];
		// need_dispatch = client_demand[t][sid_idx][c];
		if (need_dispatch == 0) continue;
		for (auto s: qos4c[c]) {
			if (! s_used[s]) continue;  // filter not used s
			t_set = &s_include_t[s];
			// dispatch to include_t = 1 ~ 5
			if ((t_set->size() < higher_check_num[s]) || t_set->find(t) != t_set->end()){  // can dispatch to 5%
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					t_set->insert(t);
					(*demand4id)[c] = 0;  // TODO: may need to delete
					break;
				}
			} else { // dispatch to < barrier(init is base_cost)  if s is used
				if (need_dispatch + ts_record[t][s] <= barrier[s]) {
					success = assign(t, s, c, sid_idx, need_dispatch);
					if (success) {
						(*demand4id)[c] = 0;  // TODO: may need to delete
						break;
					}
				}
			}
		}
	}
}

void dispatch_t_c_normal(int t, int c) {
	int need_dispatch;
	vector<int> *demand4id;
	unordered_set<int> *t_set;
	bool success;

	free_dispatch_t_c(t, c);  // free
	// pay
	for (auto &sid_idx : t_include_sid[t]) {
		demand4id = &client_demand[t][sid_idx];
		need_dispatch = (*demand4id)[c];
		int over_flow; float cost;
		if (need_dispatch == 0) continue;
		vector<float> cost_each_s;
		for (auto s: qos4c[c]){  // calculate over flow cost
			if (s_used[s]) {
				t_set = &s_include_t[s];
				if (t_set->find(t) == t_set->end()){  // in 95%, not in 5%
					int sum_here = ts_record[t][s] + need_dispatch;
					if (sum_here > bandwidth[s]) {
						cost_each_s.push_back(FLT_MAX);  // don't use it
					} else {  
						if (ts_record[t][s] < base_cost) {  // only calc the cost that over flow the base cost
							over_flow = sum_here - base_cost;
							cost_each_s.push_back(pow(over_flow, 2) / bandwidth[s] + (sum_here - base_cost));
						} else {  // already >= base_cost
							// cost = dW * ( 1 + (W'+W-2V)/c)  // W: cost at 95%, W': new cost, V: base cost, c: bandwidth
							cost = need_dispatch * (1 + (sum_here + ts_record[t][s] - 2 * base_cost) / bandwidth[s]);
							cost_each_s.push_back(cost);
						}  // WARNING: add up may < barrier, but previous step will prevent this situation
					} 
				} else {  // ???? if in 5%, will be dispatch in previous for_loop
					cost_each_s.push_back(FLT_MAX);
					success = assign(t, s, c, sid_idx, need_dispatch);
					if (success) {
						cout << "WARNING: dispatch at 5%";
						(*demand4id)[c] = 0;  // TODO: may need to delete
					}
					// TODO: try to add value
				}  
			} else {  // not used server, can as a candidate
				cost_each_s.push_back(FLT_MAX / 2);
			}
		}
		
		// open new

	}
}


void dispatch_t_c(int t, int c){
	int need_dispatch, min_idx; //, sid_idx;
	bool success;
	vector<int> *demand4id;
	unordered_set<int> *t_set;
	//////////// free to dispatch  TODO: can be a single function /////////////
	for (auto &sid_idx : t_include_sid[t]) {
		demand4id = &client_demand[t][sid_idx];
		need_dispatch = (*demand4id)[c];
		// need_dispatch = client_demand[t][sid_idx][c];
		if (need_dispatch == 0) continue;
		for (auto s: qos4c[c]) {
			if (! s_used[s]) continue;  // filter not used s
			t_set = &s_include_t[s];
			// dispatch to include_t = 1 ~ 5
			if ((t_set->size() < higher_95_num) || t_set->find(t) != t_set->end()){  // can dispatch to 5%
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					t_set->insert(t);
					(*demand4id)[c] = 0;  // TODO: may need to delete
					break;
				}
			} else { // dispatch to < barrier(init is base_cost)  if s is used
				if (need_dispatch + ts_record[t][s] <= barrier[s]) {
					success = assign(t, s, c, sid_idx, need_dispatch);
					if (success) {
						(*demand4id)[c] = 0;  // TODO: may need to delete
						break;
					}
				}
			}
		}
	}
	//////////////////////// pay to dispatch /////////////////////////////////
	int over_flow, cost;
	for (auto &sid_idx : t_include_sid[t]) {
		demand4id = &client_demand[t][sid_idx];
		need_dispatch = (*demand4id)[c];
		// need_dispatch = client_demand[t][sid_idx][c];
		if (need_dispatch == 0) continue;
		// cost of dispatch over base cost 
		vector<float> cost_each_s;
		for (auto s: qos4c[c]){  // calculate over flow cost
			if (s_used[s]) {
				t_set = &s_include_t[s];
				if (t_set->find(t) == t_set->end()){  // in 95%, not in 5%
					int sum_here = ts_record[t][s] + need_dispatch;
					if (sum_here > bandwidth[s]) {
						cost_each_s.push_back(FLT_MAX);  // don't use it
					} else {  
						if (ts_record[t][s] < base_cost) {  // only calc the cost that over flow the base cost
							over_flow = sum_here - base_cost;
							cost_each_s.push_back(pow(over_flow, 2) / bandwidth[s] + (sum_here - base_cost));
						} else {  // already >= base_cost
							// cost = dW * ( 1 + (W'+W-2V)/c)  // W: cost at 95%, W': new cost, V: base cost, c: bandwidth
							cost = need_dispatch * (1 + (sum_here + ts_record[t][s] - 2 * base_cost) / bandwidth[s]);
							cost_each_s.push_back(cost);
						}  // WARNING: add up may < barrier, but previous step will prevent this situation
					} 
				} else {  // ???? if in 5%, will be dispatch in previous for_loop
					cost_each_s.push_back(FLT_MAX);
					success = assign(t, s, c, sid_idx, need_dispatch);
					if (success) {
						cout << "WARNING: dispatch at 5%";
						(*demand4id)[c] = 0;  // TODO: may need to delete
					}
					// TODO: try to add value
				}  
			} else {  // not used server, can as a candidate
				cost_each_s.push_back(FLT_MAX / 2);
			}
		}
		min_idx = argmin(cost_each_s.begin(), cost_each_s.end());
		///////////////// over flow current server
		if (cost_each_s[min_idx] < base_cost){  
			int s = qos4c[c][min_idx];
			success = assign(t, s, c, sid_idx, need_dispatch);
			if (success) {
				(*demand4id)[c] = 0;  // TODO: may need to delete
				barrier[s] = max(barrier[s], ts_record[t][s]);
				// TODO: move other res to this server
			} else cout << "Fail at aaa";
		} else {  
		/////////////////// open new server  // TODO: check whether the bandwidth is enough
			vector<int> cand_bw, avail_num;  // choose a new s: 1. bandwidth 2. connectivity
			for (auto s: qos4c[c]){
				if (s_used[s]) cand_bw.push_back(0);
				else {
					cand_bw.push_back(bandwidth[s]);
				}
				avail_num.push_back(qos4s_num[s]);
			}
			// int idx = argmax(cand_bw.begin(), cand_bw.end());
			int idx = argmax2(cand_bw, avail_num);
			int s = qos4c[c][idx];
			if (s_used[s]) {  // s is used 
				s = qos4c[c][min_idx];  // need to overflow the server
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					(*demand4id)[c] = 0;  // TODO: may need to delete
					barrier[s] = max(barrier[s], ts_record[t][s]);
				// TODO: move other res to this server
				} else {
					cout << "Fail at bbb";
				} 
			} else {  // open a new s
				success = assign(t, s, c, sid_idx, need_dispatch);
				if (success) {
					(*demand4id)[c] = 0;  // TODO: may need to delete
					s_used[s] = true;
					t_set = &s_include_t[s];
					t_set->insert(t);
				} else cout << "Fail at ccc ";  // bandwidth of new s is not enough
				free_dispatch_a_time(t, s, c); /////////////////////////////////////////
			}
		}
	}
}

void dispatch_tc_matrix(){
	// vector<int> idxs = argsort(demand_each_t, true);
	vector<int> idxs = argpartialsort(c_demand_one_line, c_demand_one_line.size() / 3, true);
	int t, c;
	for (auto &i : idxs) {
		t = i / c_len; c = i % c_len;
		dispatch_t_c(t, c);
		// dispatch_t_c_normal(t, c);
	}
	
}

void test() { //   0   1 2 3  4 5  6 7  8 9  10 11 12
	// vector<int> a {57,22,1,67,8,34,4,58,7,89,22,72,33};
	// vector<int> idx = argpartialsort(a, 2, true);
	// int b = 1;
};

int main(int argc, char *argv[]) {
	if (argc >= 2) is_local = true;
	test();
	clock_t whole_st = clock(), st = clock();
	get_data();
	cout << "get data: \t" << (double) (clock() - st) / CLOCKS_PER_SEC << "s" << endl;

	st = clock();
	analyse_qos();
	cout << "analyse_qos: \t" << (double) (clock() - st) / CLOCKS_PER_SEC << "s" << endl;

	st = clock();
	// dispatch_baseline();
	// dispatch_top5_opt();
	dispatch_tc_matrix();
	cout << "dispatch: \t" << (double) (clock() - st) / CLOCKS_PER_SEC << "s" << endl;

	st = clock();
	output_result();
	cout << "output: \t" << (double) (clock() - st) / CLOCKS_PER_SEC << "s" << endl;

	st = clock();
	cout << "whole time: \t" << (double) (clock() - whole_st) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}

