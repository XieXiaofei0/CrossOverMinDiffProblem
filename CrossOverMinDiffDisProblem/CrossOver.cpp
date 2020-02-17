#include"CrossOver.h"
#include"utility.h"
#include"InputOutput.h"
#include <queue>

using namespace std;
using namespace xxf_utility;

namespace min_diff_dp {

    LogSwitch logsw_local(1, 1, "PINCOL_PT");

    bool compareByAscend(const pair<int, Distance> &a, const  pair<int, Distance> &b) {
        return a.second < b.second;
    }
    bool compareByDescend(const pair<int, Distance> &a, const  pair<int, Distance> &b) {
        return a.second > b.second;
    }

    CrossOver::CrossOver(const UMatrix &_matrix, int _tabu_iter, int _tabu_length, int popu, double param, const Solution &_init_sol_one, const Solution &_init_sol_two,
       const int _size_of_tabu, const double _param1, const double _param2, const double _param3) :
        ins(_matrix), nb_nodes(_matrix.setele_num()), nb_sub_nodes(_matrix.subsetele_num()),
        population(popu), single_iter_length(_tabu_iter), tabu_iter_max(_tabu_length), generation(0),
        max_time(nb_nodes * 1.5 * 1000), rate_of_sele_nodes(param),
        size_of_tabu_list(_size_of_tabu), hashFun_one_param(_param1), hashFun_two_param(_param2), hashFun_three_param(_param3)
    {
        tabu_list_one.resize(size_of_tabu_list, 0);
        tabu_list_two.resize(size_of_tabu_list, 0);
        tabu_list_three.resize(size_of_tabu_list, 0);
        //Ҳ���Լ��������
        cross_best.resize(nb_nodes);
        node_value.resize(population, vector<int>(nb_nodes));
        node_value[0] = _init_sol_one.node_values();
        node_value[1] = _init_sol_two.node_values();
        cur_obj.resize(2);
        local_best.resize(population, vector<int>(nb_nodes));
        local_best[0] = _init_sol_one.node_values();
        local_best[1] = _init_sol_two.node_values();
        local_best_obj.resize(2);
        hash_key_temp_one.resize(nb_nodes, 0);
        hash_key_temp_two.resize(nb_nodes, 0);
        hash_key_temp_three.resize(nb_nodes, 0);
        init();
    }

    void CrossOver::init() {
        size_neighbor_struc = (int)((nb_nodes - nb_sub_nodes)*rate_of_sele_nodes);   //��ʼ��List��С�Ͳ���

        node_dis_sum.resize(population);
        best_solu_dis_sum.resize(population);
        select_nodes.resize(population);
        no_select_nodes.resize(population);
        max_select_node.resize(2);
        min_select_node.resize(2);
        best_max_select_node.resize(2);
        best_min_select_node.resize(2);
        best_max_select_node[0] = -1;
        best_max_select_node[1] = -1;
        best_min_select_node[0] = DISTANCE_MAX;
        best_min_select_node[1] = DISTANCE_MAX;
        best_hashfun_one.resize(2);
        best_hashfun_two.resize(2);
        best_hashfun_three.resize(2);
        for (int i = 0; i < population; ++i) {
            node_dis_sum[i].reserve(nb_nodes);
            best_solu_dis_sum[i].reserve(nb_nodes);
            select_nodes[i].reserve(nb_sub_nodes);
            no_select_nodes[i].reserve(nb_nodes - nb_sub_nodes);
        }
        List<long long> sum1(2, 0), sum2(2, 0), sum3(2, 0);
        for (int i = 0; i < nb_nodes; ++i) {
            int one = (int)(floor(pow(i, hashFun_one_param)));
            int two = (int)(floor(pow(i, hashFun_two_param)));
            int three = (int)(floor(pow(i, hashFun_three_param)));
            hash_key_temp_one[i] = one;
            hash_key_temp_two[i] = two;
            hash_key_temp_three[i] = three;
            for (int j = 0; j < population; ++j) {
                double dis = 0.0, dis_elite = 0.0;
                for (int k = 0; k < nb_nodes; ++k) {
                    if (node_value[j][k])dis += ins.dis_nodes(i, k);
                }
                node_dis_sum[j].push_back(dis);
                best_solu_dis_sum[j].push_back(dis);
                if (node_value[j][i]) {
                    if (max_select_node[j] < dis)max_select_node[j] = dis;
                    if (min_select_node[j] > dis)min_select_node[j] = dis;
                    sum1[j] += hash_key_temp_one[i];                      //���㵱ǰ��ʷ���Ž�������ϣ�����ļ�ֵ
                    sum2[j] += hash_key_temp_two[i];
                    sum3[j] += hash_key_temp_three[i];
                }
            }
        }
        best_hashfun_one[0] = sum1[0] % size_of_tabu_list;
        best_hashfun_one[1] = sum1[1] % size_of_tabu_list;
        best_hashfun_two[0] = sum2[0] % size_of_tabu_list;
        best_hashfun_two[1] = sum2[1] % size_of_tabu_list;
        best_hashfun_three[0] = sum3[0] % size_of_tabu_list;
        best_hashfun_three[1] = sum3[1] % size_of_tabu_list;
        cur_obj[0] = max_select_node[0] - min_select_node[0];
        cur_obj[1] = max_select_node[1] - min_select_node[1];
        local_best_obj[0] = cur_obj[0];
        local_best_obj[1] = cur_obj[1];

        if (local_best_obj[0] < local_best_obj[1]) {
            cross_best = local_best[0];
            cross_best_obj = local_best_obj[0];
        }
        else {
            cross_best = local_best[1];
            cross_best_obj = local_best_obj[1];
        }

        for (int i = 0; i < population; ++i) {
            Distance temp_sum = (max_select_node[i] + min_select_node[i]) / 2;
            for (int j = 0; j < nb_nodes; ++j) {                                     //��ʼ�����������ṹno_select_nodes��select_nodes
                Distance temp = fabs(node_dis_sum[i][j] - temp_sum);
                if (node_value[i][j])select_nodes[i].push_back(make_pair(j, temp));
                else no_select_nodes[i].push_back(make_pair(j, temp));
            }
            sort(no_select_nodes[i].begin(), no_select_nodes[i].end(), compareByAscend);
        }
    }

    Solution CrossOver::solve() {
        Timer time(max_time);
        while (!time.isTimeOut()) {
            //ͨ��cross�е��жϣ�����������ƶȲ���Խ��Խ��
            //��������н���õ�������
            cross();
            //������ֱ���н�������
            tabusearch();
            //�������Ž�
            update();
            generation++;
        }
        return Solution(nb_nodes, nb_sub_nodes, cross_best, cross_best_obj);
    }

    void CrossOver::cross() {
        List<List<int>> solu;
        solu.resize(2);
        solu[0].reserve(nb_sub_nodes);
        solu[1].reserve(nb_sub_nodes);
        List<List<int>> cur_nodes(2, List<int>(nb_nodes, 0));
        List<int> hashone(2, 0), hashtwo(2, 0), hashthree(2, 0);
        for (int i = 0; i < population; ++i) {
            for (int j = 0; j < nb_nodes; ++j) {
                if (local_best[i][j])solu[i].push_back(j);
            }
        }
        //��������ڵ�
        sort(solu[0].begin(), solu[0].end());
        sort(solu[1].begin(), solu[1].end());
        vector<int> ivec(nb_sub_nodes);
        auto iter = set_intersection(solu[0].begin(), solu[0].end(), solu[1].begin(), solu[1].end(), ivec.begin());
        ivec.resize(iter - ivec.begin());//����ȷ��ivec��С
        for (int i = 0; i < ivec.size(); ++i) {
            cur_nodes[0][ivec[i]] = 1;
            cur_nodes[1][ivec[i]] = 1;
            hashone[0] += hash_key_temp_one[ivec[i]];      //ͬʱ�����µĽ��hashֵ
            hashtwo[0] += hash_key_temp_two[ivec[i]];
            hashthree[0] += hash_key_temp_three[ivec[i]];
            hashone[1] = hashone[0];
            hashtwo[1] = hashtwo[0];
            hashthree[1] = hashthree[0];
        }

        //test
        //mylog << "������������Ž⣺" <<= logsw_info;
        //mylog << "��1��  " << local_best_obj[0] << "     ��2��  " << local_best_obj[1] <<= logsw_info;
        //mylog << "��ǰ�ǵڼ�����" << cycle << "      ��ǰ�������Ϊ��" << generation << "     ��ͬ�ڵ������" << ivec.size() <<= logsw_info;
        //test end

        vector<vector<int>> select_nodes;
        select_nodes.resize(2);
        select_nodes[0] = ivec;
        select_nodes[0].reserve(nb_sub_nodes);
        select_nodes[1] = ivec;
        select_nodes[1].reserve(nb_sub_nodes);
        int remain = nb_sub_nodes - ivec.size();
        //�����ѡʣ��Ľڵ�
        srand(myrand.getSeed());
        int cur = 0;
        for (int i = 0; i < population; ++i) {
            int num = 0, index = 0;
            while (num < remain) {
                while (true) {
                    index = rand() % nb_sub_nodes;
                    if (cur_nodes[i][solu[cur][index]] == 0)break;
                }
                cur_nodes[i][solu[cur][index]] = 1;
                select_nodes[i].push_back(solu[cur][index]);
                hashone[i] += hash_key_temp_one[solu[cur][index]];
                hashtwo[i] += hash_key_temp_two[solu[cur][index]];
                hashthree[i] += hash_key_temp_three[solu[cur][index]];
                cur = 1 - cur;
                ++num;
            }
            best_hashfun_one[i] = hashone[i] % size_of_tabu_list;
            best_hashfun_two[i] = hashtwo[i] % size_of_tabu_list;
            best_hashfun_three[i] = hashthree[i] % size_of_tabu_list;
        }
        local_best[0] = cur_nodes[0];
        local_best[1] = cur_nodes[1];

        //test
        //mylog << "��һ����Ϊ��  " <<= logsw_info;
        //for (int i = 0; i < nb_nodes; ++i)if (local_best[0][i])mylog << i << " ";
        //mylog <<= logsw_info;
        //mylog << "�ڶ�����Ϊ��   " <<= logsw_info;
        //for (int i = 0; i < nb_nodes; ++i)if (local_best[1][i])mylog << i << " ";
        //sort(select_nodes[0].begin(), select_nodes[0].end());
        //sort(select_nodes[1].begin(), select_nodes[1].end());
        //vector<int> ivec1(nb_sub_nodes);
        //auto iter1 = set_intersection(select_nodes[0].begin(), select_nodes[0].end(), select_nodes[1].begin(), select_nodes[1].end(), ivec1.begin());
        //ivec1.resize(iter1 - ivec1.begin());//����ȷ��ivec��С

        //test end

        for (int i = 0; i < population; ++i) {
            Distance max = -1, min = DISTANCE_MAX;
            for (int l = 0; l < nb_nodes; ++l) {
                best_solu_dis_sum[i][l] = 0;
                for (int m = 0; m < nb_sub_nodes; ++m) {
                    best_solu_dis_sum[i][l] += ins.dis_nodes(l, select_nodes[i][m]);
                }
                if (local_best[i][l]) {
                    if (max < best_solu_dis_sum[i][l])max = best_solu_dis_sum[i][l];
                    if (min > best_solu_dis_sum[i][l])min = best_solu_dis_sum[i][l];
                }
            }
            best_max_select_node[i] = max;
            best_min_select_node[i] = min;
            local_best_obj[i] = max - min;
        }
        //test
        //mylog << "������֮��Ľ�Ϊ�� " <<= logsw_info;
        //mylog << "��1��  " << local_best_obj[0] << "     ��2��  " << local_best_obj[1] <<= logsw_info;
        //mylog << "��ͬ�ڵ����Ϊ��    " << ivec1.size() <<= logsw_info;
        //mylog <<= logsw_info;
        //mylog <<= logsw_info;
        //test end
        //test
        //List<int> select;
        //select.reserve(nb_sub_nodes);
        //for (int i = 0; i < nb_nodes; ++i) {
        //    if (local_best[0][i])select.push_back(i);
        //}
        //Distance min = DISTANCE_MAX, max = 0;
        //for (int i = 0; i < nb_sub_nodes; ++i) {
        //    Distance cur = 0;
        //    for (int j = 0; j < nb_sub_nodes; ++j) {
        //        cur += ins.dis_nodes(select[i], select[j]);
        //    }
        //    if (cur > max)max = cur;
        //    if (cur < min)min = cur;
        //}
        //mylog << "���ֵ��" << max << "     " << best_max_select_node[0] <<= logsw_info;
        //Distance obj = max - min;
        //mylog << "����Ŀ�꺯��ֵΪ��" << obj <<= logsw_info;
        //mylog << "��ǰ��Ŀ�꺯��ֵΪ��" << local_best_obj[0] <<= logsw_info;
        //long long sum = 0;
        //for (int i = 0; i < nb_nodes; ++i) {
        //    if (local_best[1][i]) {
        //        sum += (int)(floor(pow(i, hashFun_one_param)));
        //    }
        //}
        //mylog << "��ϣֵ��  " << sum << "    " << best_hashfun_one[1] <<= logsw_info;
        //test end

    }

    void CrossOver::tabusearch() {
        //TODO�����Զ��ٲ�֮�ڸĽ�������ʷ���Ž�ֹͣ��Ŀǰ����Ϊһ������
        for (int i = 0; i < population; ++i) {
            bool tabu_flag = false;
            int iter = 0;
            int step_length = 0;         //������ٲ�֮�ڸĽ�������ʷ���Ž�
            //while (iter < tabu_iter_max) {
            while(true){
                int _hashfun_one = best_hashfun_one[i];
                int _hashfun_two = best_hashfun_two[i];
                int _hashfun_three = best_hashfun_three[i];
                node_value[i] = local_best[i];               //ǿ���������ԣ�����ʷ���Ž���µ�ǰ��
                cur_obj[i] = local_best_obj[i];
                //����ʷ���Ž�!=��ǰ��,��Ҫ���µ�ǰ��Ĳ������ݽṹnode_dis_sum��max_select_node�� min_select_node��select_nodes��no_select_nodes
                max_select_node[i] = best_max_select_node[i];                   //xxf:�����bug1��֮ǰ�滻��ǰ��ʱ��δ����������ݽṹnode_dis_sum��select_nodes��no_select_nodes
                min_select_node[i] = best_min_select_node[i];
                no_select_nodes[i].clear();            //���¸����ṹselect_nodes��no_select_nodes
                select_nodes[i].clear();
                Distance temp_sum = (best_max_select_node[i] + best_min_select_node[i]) / 2;
                for (int j = 0; j < nb_nodes; ++j) {
                    node_dis_sum[i][j] = best_solu_dis_sum[i][j];
                    Distance temp = fabs(best_solu_dis_sum[i][j] - temp_sum);
                    if (node_value[i][j])select_nodes[i].push_back(make_pair(j, temp));
                    else no_select_nodes[i].push_back(make_pair(j, temp));
                }
                sort(no_select_nodes[i].begin(), no_select_nodes[i].end(), compareByAscend);    //δѡ�е���������

                int count = 0;

                while (count <= single_iter_length) {
                    pair<int, int> swap_pair(-1, -1);              //����ṹ������õķǽ��ɽ�:����ǽ��ɵ���õĽ�����;��һ��I1-->I0���ڶ���I0-->I1
                    pair<Distance, Distance> new_obj(DISTANCE_MAX, 0);         //�����Ӧ��Ŀ�꺯�������������С����
                    find_best_move(i, swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three);   //_hashfun_one�Ѿ��Ǹ��½�Ĺ�ϣ����ֵ
                    //TODO:�ж��Ƿ���������ⶼ�ڽ�����
                    if (swap_pair.first == -1) {       //�ж��Ƿ�����ⶼ�ڽ�����
                        mylog << "����ⶼ�ڽ����У���ǰiter��" << iter <<= logsw_local;
                        tabu_flag = true;
                        break;
                    }
                    if (update_solu(i, swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three, step_length))count = 0;   //���µ�ǰ�⡢��ʷ���Ž⡢count
                    else count++;
                    iter++;
                    if (step_length == tabu_iter_max) {
                        //iter_max += 10;
                        tabu_flag = true;
                        break;
                    }
                }
                if (tabu_flag)break;
            }

        }
    }

    void CrossOver::update() {
        int min_index = -1;
        if (local_best_obj[0] < local_best_obj[1]) {
            min_index = 0;
        }
        else min_index = 1;
        if (local_best_obj[min_index] < cross_best_obj) {
            cross_best = local_best[min_index];
            cross_best_obj = local_best_obj[min_index];
            mylog << "    ��ǰ�������Ϊ:    " << generation << "       �Ľ���ʷ���Ž⣺  " << cross_best_obj <<= logsw_info;
        }
    }

    void CrossOver::find_best_move(const int number, pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three) {   //xxf:done,right--12.10
        int new_hashone = _hash_one;     //���浱ǰ��Ĺ�ϣ����ֵ
        int new_hashtwo = _hash_two;
        int new_hashthree = _hash_three;
        for (int i = 0; i < nb_sub_nodes; ++i)
        {
            int one_toZero_node = select_nodes[number][i].first;
            for (int j = 0; j < size_neighbor_struc; ++j) {
                int zero_toOne_node = no_select_nodes[number][j].first;
                int _new_hash_one = new_hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //���㽻�����½�Ĺ�ϣ����ֵ
                int _new_hash_two = new_hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
                int _new_hash_three = new_hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
                _new_hash_one = (_new_hash_one + size_of_tabu_list) % size_of_tabu_list;                //xxf�����bug2����ֹ���ָ�����>size_of_tabu_list����
                _new_hash_two = (_new_hash_two + size_of_tabu_list) % size_of_tabu_list;
                _new_hash_three = (_new_hash_three + size_of_tabu_list) % size_of_tabu_list;
                if (tabu_list_three[_new_hash_three]) {       //xxf��������������ϣ����Ϊ�˼�С��ͻ��ײ��һ����ϣ�������ײ�ͬ��ӳ�䵽ͬһ��valueֵ���Һ����׳�������ⶼ�������������Ϊ�Ǵ�ǿ���������Կ�ʼ
                    if (tabu_list_two[_new_hash_two])
                        if (tabu_list_one[_new_hash_one]) {
                            continue;
                        }
                }
                Distance temp_min = node_dis_sum[number][zero_toOne_node] - ins.dis_nodes(zero_toOne_node, one_toZero_node);      //��¼��ǰһ�ζ��������ֵ����Сֵ;��ʼ��Ϊ����֮���µ�ѡ�нڵ�ľ���֮��
                Distance temp_max = temp_min;
                for (int k = 0; k < nb_sub_nodes; ++k) {
                    if (k == i)continue;
                    int node = select_nodes[number][k].first;
                    Distance update_dis = node_dis_sum[number][node] - ins.dis_nodes(node, one_toZero_node) + ins.dis_nodes(node, zero_toOne_node);
                    if (update_dis > temp_max)temp_max = update_dis;
                    else if (temp_min > update_dis)temp_min = update_dis;
                    else;
                }
                if ((_new_obj.first - _new_obj.second) > (temp_max - temp_min)) {   //����������
                    _pair.first = one_toZero_node;
                    _pair.second = zero_toOne_node;
                    _new_obj.first = temp_max;
                    _new_obj.second = temp_min;
                    _hash_one = _new_hash_one;
                    _hash_two = _new_hash_two;
                    _hash_three = _new_hash_three;
                }
            }
        }
    }

    bool CrossOver::update_solu(const int number, const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, int &step) {  //xxf:done,right-12.10
        bool flag = false;         //��ʾ�Ƿ������ʷ���Ž�
        node_value[number][_pair.first] = 0;                       //���µ�ǰ��
        node_value[number][_pair.second] = 1;
        cur_obj[number] = _new_obj.first - _new_obj.second;
        max_select_node[number] = _new_obj.first;
        min_select_node[number] = _new_obj.second;
        if (local_best_obj[number] > cur_obj[number])       //������ʷ���Ž⣬������ʷ���Ž����ؽṹ
        {
            local_best[number] = node_value[number];  //����ܸĽ���ʷ���Ž⣬�������ʷ���Ž�,����true
            local_best_obj[number] = cur_obj[number];
            best_hashfun_one[number] = _hash_one;
            best_hashfun_two[number] = _hash_two;
            best_hashfun_three[number] = _hash_three;
            //test:TODO:���ٲ�֮���ܵ�������
            //mylog << "��" << number << " ��popu  ��\n��ǰΪ��" << local_best_obj[number] << "     ������� ��" << cross_num << "    ��������  �� " <<= logsw_local;
            //test end
            flag = true;
            step = 0;
        }
        step++;
        tabu_list_one[_hash_one] = 1;         //�������������б�
        tabu_list_two[_hash_two] = 1;
        tabu_list_three[_hash_three] = 1;
        no_select_nodes[number].clear();           //���¸����ṹselect_nodes��no_select_nodes
        select_nodes[number].clear();
        Distance temp_sum = (max_select_node[number] + min_select_node[number]) / 2.0;
        if (flag) {           //����������ʷ���Ž⣬�򱣴���ʷ���Ž��������ݽṹbest_solu_dis_sum
            best_max_select_node[number] = max_select_node[number];
            best_min_select_node[number] = min_select_node[number];
            for (int i = 0; i < nb_nodes; ++i) {
                node_dis_sum[number][i] = node_dis_sum[number][i] - ins.dis_nodes(i, _pair.first) + ins.dis_nodes(i, _pair.second);
                best_solu_dis_sum[number][i] = node_dis_sum[number][i];
                Distance temp = fabs(node_dis_sum[number][i] - temp_sum);
                if (node_value[number][i])select_nodes[number].push_back(make_pair(i, temp));
                else no_select_nodes[number].push_back(make_pair(i, temp));
            }
        }
        else {
            for (int i = 0; i < nb_nodes; ++i) {
                node_dis_sum[number][i] = node_dis_sum[number][i] - ins.dis_nodes(i, _pair.first) + ins.dis_nodes(i, _pair.second);
                Distance temp = fabs(node_dis_sum[number][i] - temp_sum);
                if (node_value[number][i])select_nodes[number].push_back(make_pair(i, temp));
                else no_select_nodes[number].push_back(make_pair(i, temp));
            }
        }
        sort(no_select_nodes[number].begin(), no_select_nodes[number].end(), compareByAscend);    //δѡ�е���������
        return flag;
    }
}