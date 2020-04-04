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

    CrossOver::CrossOver(const UMatrix &_matrix, int _tabu_iter, int _tabu_length, int _iter_cycle, int popu, double param, const Solution &_init_sol_one, const Solution &_init_sol_two,
        const Solution &_init_elite_one, const Solution &_init_elite_two, const int _size_of_tabu, const double _param1, const double _param2, const double _param3) :
        ins(_matrix), iter_cycle(_iter_cycle), cycle(0), nb_nodes(_matrix.setele_num()), nb_sub_nodes(_matrix.subsetele_num()),
        population(popu), single_iter_length(_tabu_iter), tabu_iter_max(_tabu_length), generation(0),
        max_time(nb_nodes * 1.5 * 1000), rate_of_sele_nodes(param),
        size_of_tabu_list(_size_of_tabu), hashFun_one_param(_param1), hashFun_two_param(_param2), hashFun_three_param(_param3)
    {
        tabu_list_one.resize(size_of_tabu_list, 0);
        tabu_list_two.resize(size_of_tabu_list, 0);
        tabu_list_three.resize(size_of_tabu_list, 0);
        //也可自己随机生成
        cross_best.resize(nb_nodes);
        elite.resize(population, vector<int>(nb_nodes));
        elite[0] = _init_elite_one.node_values();
        elite[1] = _init_elite_two.node_values();
        elite_obj.resize(2);
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
        size_neighbor_struc = (int)((nb_nodes - nb_sub_nodes)*rate_of_sele_nodes);   //初始化List大小和参数
        //临时数据结构
        List<List<Distance>> elite_dis_sum;
        elite_dis_sum.resize(population);
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
            elite_dis_sum[i].reserve(nb_nodes);
            node_dis_sum[i].reserve(nb_nodes);
            best_solu_dis_sum[i].reserve(nb_nodes);
            select_nodes[i].reserve(nb_sub_nodes);
            no_select_nodes[i].reserve(nb_nodes - nb_sub_nodes);
        }
        List<long long> sum1(2, 0), sum2(2, 0), sum3(2, 0);
        //临时数据结构：记录elite的目标函数值
        List<Distance> elite_max_select_node(2, -1);
        List<Distance> elite_min_select_node(2, DISTANCE_MAX);
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
                    if (elite[j][k])dis_elite += ins.dis_nodes(i, k);
                }
                node_dis_sum[j].push_back(dis);
                best_solu_dis_sum[j].push_back(dis);
                elite[j].push_back(dis_elite);
                if (node_value[j][i]) {
                    if (max_select_node[j] < dis)max_select_node[j] = dis;
                    if (min_select_node[j] > dis)min_select_node[j] = dis;
                    sum1[j] += hash_key_temp_one[i];                      //计算当前历史最优解三个哈希函数的键值
                    sum2[j] += hash_key_temp_two[i];
                    sum3[j] += hash_key_temp_three[i];
                }
                if (elite[j][i]) {
                    if (elite_max_select_node[j] < dis_elite)elite_max_select_node[j] = dis_elite;
                    if (elite_min_select_node[j] > dis_elite)elite_min_select_node[j] = dis_elite;
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
        elite_obj[0] = elite_max_select_node[0] - elite_min_select_node[0];
        elite_obj[1] = elite_max_select_node[1] - elite_min_select_node[1];
        local_best_obj[0] = cur_obj[0];
        local_best_obj[1] = cur_obj[1];

        if (cur_obj[0] < cur_obj[1]) {
            if (elite_obj[0] < elite_obj[1]) {
                if (cur_obj[0] < elite_obj[0]) {
                    cross_best = node_value[0];
                    cross_best_obj = cur_obj[0];
                }
                else {
                    cross_best = elite[0];
                    cross_best_obj = elite_obj[0];
                }
            }
            else {
                if (cur_obj[0] < elite_obj[1]) {
                    cross_best = node_value[0];
                    cross_best_obj = cur_obj[0];
                }
                else {
                    cross_best = elite[1];
                    cross_best_obj = elite_obj[1];
                }
            }
        }
        else {
            if (elite_obj[0] < elite_obj[1]) {
                if (cur_obj[1] < elite_obj[0]) {
                    cross_best = node_value[1];
                    cross_best_obj = cur_obj[1];
                }
                else {
                    cross_best = elite[0];
                    cross_best_obj = elite_obj[0];
                }
            }
            else {
                if (cur_obj[1] < elite_obj[1]) {
                    cross_best = node_value[1];
                    cross_best_obj = cur_obj[1];
                }
                else {
                    cross_best = elite[1];
                    cross_best_obj = elite_obj[1];
                }
            }
        }
        for (int i = 0; i < population; ++i) {
            Distance temp_sum = (max_select_node[i] + min_select_node[i]) / 2;
            for (int j = 0; j < nb_nodes; ++j) {                                     //初始化两个辅助结构no_select_nodes和select_nodes
                Distance temp = fabs(node_dis_sum[i][j] - temp_sum);
                if (node_value[i][j])select_nodes[i].push_back(make_pair(j, temp));
                else no_select_nodes[i].push_back(make_pair(j, temp));
            }
            sort(no_select_nodes[i].begin(), no_select_nodes[i].end(), compareByAscend);
        }
    }

    Solution CrossOver::solve() {
        Timer time(max_time);
        clock_t start_time = clock();
        while (!time.isTimeOut()) {
            //通过cross中的判断，两个解的相似度不会越来越大
            //两个解进行交叉得到两个解
            cross();
            //两个解分别进行禁忌搜索
            tabusearch();
            //更新最优解和elite1解
            update_elite();
            update_best();
            if ((generation%iter_cycle) == 0) {
                local_best[0] = elite[1];
                local_best_obj[0] = elite_obj[1];
                elite[1] = elite[0];
                elite_obj[1] = elite_obj[0];
                random_elite_one();
                cycle++;
            }
            generation++;
        }
        mylog << "\n总交叉代数：" << generation <<= logsw_local;
        clock_t end_time = clock();
        //test
        //关键信息输出到文件中
        ofstream outFile;
        outFile.open("../Deploy/Logs/log.csv", ios::app);
        outFile << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << ',' << nb_nodes << ',' << nb_sub_nodes << ',' << myrand.getSeed() << ',' << tabu_iter_max << ',' << iter_cycle << ',' << generation << ',' << cross_best_obj << endl;
        outFile.close();
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
        //求出公共节点
        sort(solu[0].begin(), solu[0].end());
        sort(solu[1].begin(), solu[1].end());
        vector<int> ivec(nb_sub_nodes);
        auto iter = set_intersection(solu[0].begin(), solu[0].end(), solu[1].begin(), solu[1].end(), ivec.begin());
        ivec.resize(iter - ivec.begin());//重新确定ivec大小
        for (int i = 0; i < ivec.size(); ++i) {
            cur_nodes[0][ivec[i]] = 1;
            cur_nodes[1][ivec[i]] = 1;
            hashone[0] += hash_key_temp_one[ivec[i]];      //同时保存新的解的hash值
            hashtwo[0] += hash_key_temp_two[ivec[i]];
            hashthree[0] += hash_key_temp_three[ivec[i]];
            hashone[1] = hashone[0];
            hashtwo[1] = hashtwo[0];
            hashthree[1] = hashthree[0];
        }
        vector<vector<int>> select_nodes;
        select_nodes.resize(2);
        select_nodes[0] = ivec;
        select_nodes[0].reserve(nb_sub_nodes);
        select_nodes[1] = ivec;
        select_nodes[1].reserve(nb_sub_nodes);
        int remain = nb_sub_nodes - ivec.size();
        //随机挑选剩余的节点
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
    }

    void CrossOver::cross_over() {
        //求出两个解的公共节点
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
        sort(solu[0].begin(), solu[0].end());
        sort(solu[1].begin(), solu[1].end());
        vector<int> ivec(nb_sub_nodes);
        auto iter = set_intersection(solu[0].begin(), solu[0].end(), solu[1].begin(), solu[1].end(), ivec.begin());
        ivec.resize(iter - ivec.begin());//重新确定ivec大小
        //end
        //找出两个解的非公共节点放在两个数组ivec_difference_one和ivec_difference_two中
        vector<int> ivec_difference_one(nb_sub_nodes);
        auto iter_difference_one = set_difference(solu[0].begin(), solu[0].end(), ivec.begin(), ivec.end(), ivec_difference_one.begin());
        ivec_difference_one.resize(iter_difference_one - ivec_difference_one.begin());//重新确定ivec大小
        vector<int> ivec_difference_two(nb_sub_nodes);
        auto iter_difference_two = set_difference(solu[1].begin(), solu[1].end(), ivec.begin(), ivec.end(), ivec_difference_two.begin());
        ivec_difference_two.resize(iter_difference_two - ivec_difference_two.begin());//重新确定ivec大小
        //set<int> no_select_nodes_one_first(ivec_difference_one.begin(), ivec_difference_one.end());
        //set<int> no_select_nodes_one_second(ivec_difference_two.begin(), ivec_difference_two.end());
        //set<int> no_select_nodes_two_first(ivec_difference_one.begin(), ivec_difference_one.end());
        //set<int> no_select_nodes_two_second(ivec_difference_two.begin(), ivec_difference_two.end());
        //公共节点的相关结构更新；将选中节点放在集合select_nodes_one中
        set<int> select_nodes_one;
        set<int> select_nodes_two;
        for (int l = 0; l < nb_nodes; ++l) {
            best_solu_dis_sum[0][l] = 0;
            best_solu_dis_sum[1][l] = 0;
            for (int i = 0; i < ivec.size(); ++i) {
                best_solu_dis_sum[0][l] += ins.dis_nodes(l, ivec[i]);
                best_solu_dis_sum[1][l] += ins.dis_nodes(l, ivec[i]);
                cur_nodes[0][ivec[i]] = 1;
                cur_nodes[1][ivec[i]] = 1;
                select_nodes_one.insert(ivec[i]); //new
                select_nodes_two.insert(ivec[i]); //new
                hashone[0] += hash_key_temp_one[ivec[i]];      //同时保存新的解的hash值
                hashtwo[0] += hash_key_temp_two[ivec[i]];
                hashthree[0] += hash_key_temp_three[ivec[i]];
                hashone[1] = hashone[0];
                hashtwo[1] = hashtwo[0];
                hashthree[1] = hashthree[0];
            }
        }
        //end
        int remain = nb_sub_nodes - ivec.size();
        //首先交叉得出第一个解
        int index = 0;
        int num = 0;
        while (num < remain) {
            pair<int, Distance> best_node(-1, DISTANCE_MAX);   //保存最好的将要选中的节点
            Distance min_node = DISTANCE_MAX;
            Distance max_dis = -1;
            Distance min_dis = DISTANCE_MAX;
            Distance dis = 0;
            if (index == 0) {
                for (int j = 0; j < remain; ++j) {
                    if (cur_nodes[0][ivec_difference_one[j]])continue;
                    else {
                        set<int>::iterator it = select_nodes_one.begin();
                        for (; it != select_nodes_one.end(); it++) {
                            int distance = ins.dis_nodes(ivec_difference_one[j], *it);
                            dis += distance;
                            if (max_dis < best_solu_dis_sum[0][*it] + distance)max_dis = best_solu_dis_sum[0][*it] + distance;
                            if (min_dis > best_solu_dis_sum[0][*it] + distance)min_dis = best_solu_dis_sum[0][*it] + distance;
                        }
                        max_dis = max(max_dis, dis);
                        min_dis = min(min_dis, dis);
                        if (max_dis - min_dis < min_node) {
                            min_node = max_dis - min_dis;
                            best_node.first = ivec_difference_one[j];
                            best_node.second = dis;
                        }
                    }
                }
            }
            else {
                for (int j = 0; j < remain; ++j) {
                    if (cur_nodes[0][ivec_difference_two[j]])continue;
                    else {
                        set<int>::iterator it = select_nodes_one.begin();
                        for (; it != select_nodes_one.end(); it++) {
                            int distance = ins.dis_nodes(ivec_difference_two[j], *it);
                            dis += distance;
                            if (max_dis < best_solu_dis_sum[0][*it] + distance)max_dis = best_solu_dis_sum[0][*it] + distance;
                            if (min_dis > best_solu_dis_sum[0][*it] + distance)min_dis = best_solu_dis_sum[0][*it] + distance;
                        }
                        max_dis = max(max_dis, dis);
                        min_dis = min(min_dis, dis);
                        if (max_dis - min_dis < min_node) {
                            min_node = max_dis - min_dis;
                            best_node.first = ivec_difference_two[j];
                            best_node.second = dis;
                        }
                    }
                }

            }
            index = 1 - index;
            ++num;
            select_nodes_one.insert(best_node.first);
            cur_nodes[0][best_node.first] = 1;
            hashone[0] += hash_key_temp_one[best_node.first];
            hashtwo[0] += hash_key_temp_two[best_node.first];
            hashthree[0] += hash_key_temp_three[best_node.first];
            Distance max = -1, min = DISTANCE_MAX;
            for (int l = 0; l < nb_nodes; ++l) {
                best_solu_dis_sum[0][l] += ins.dis_nodes(l, best_node.first);
                if (cur_nodes[0][l]) {
                    if (max < best_solu_dis_sum[0][l])max = best_solu_dis_sum[0][l];
                    if (min > best_solu_dis_sum[0][l])min = best_solu_dis_sum[0][l];
                }
            }
            best_max_select_node[0] = max;
            best_min_select_node[0] = min;
            local_best_obj[0] = max - min;
        }
        best_hashfun_one[0] = hashone[0] % size_of_tabu_list;
        best_hashfun_two[0] = hashtwo[0] % size_of_tabu_list;
        best_hashfun_three[0] = hashthree[0] % size_of_tabu_list;
        local_best[0] = cur_nodes[0];
        //交叉得出第2个解
        index = 1;
        num = 0;
        while (num < remain) {
            pair<int, Distance> best_node(-1, DISTANCE_MAX);   //保存最好的将要选中的节点
            Distance min_node = DISTANCE_MAX;
            Distance max_dis = -1;
            Distance min_dis = DISTANCE_MAX;
            Distance dis = 0;
            if (index == 0) {
                for (int j = 0; j < remain; ++j) {
                    if (cur_nodes[1][ivec_difference_one[j]])continue;
                    else {
                        set<int>::iterator it = select_nodes_two.begin();
                        for (; it != select_nodes_two.end(); it++) {
                            int distance = ins.dis_nodes(ivec_difference_one[j], *it);
                            dis += distance;
                            if (max_dis < best_solu_dis_sum[1][*it] + distance)max_dis = best_solu_dis_sum[1][*it] + distance;
                            if (min_dis > best_solu_dis_sum[1][*it] + distance)min_dis = best_solu_dis_sum[1][*it] + distance;
                        }
                        max_dis = max(max_dis, dis);
                        min_dis = min(min_dis, dis);
                        if (max_dis - min_dis < min_node) {
                            min_node = max_dis - min_dis;
                            best_node.first = ivec_difference_one[j];
                            best_node.second = dis;
                        }
                    }
                }
            }
            else {
                for (int j = 0; j < remain; ++j) {
                    if (cur_nodes[1][ivec_difference_two[j]])continue;
                    else {
                        set<int>::iterator it = select_nodes_two.begin();
                        for (; it != select_nodes_two.end(); it++) {
                            int distance = ins.dis_nodes(ivec_difference_two[j], *it);
                            dis += distance;
                            if (max_dis < best_solu_dis_sum[1][*it] + distance)max_dis = best_solu_dis_sum[1][*it] + distance;
                            if (min_dis > best_solu_dis_sum[1][*it] + distance)min_dis = best_solu_dis_sum[1][*it] + distance;
                        }
                        max_dis = max(max_dis, dis);
                        min_dis = min(min_dis, dis);
                        if (max_dis - min_dis < min_node) {
                            min_node = max_dis - min_dis;
                            best_node.first = ivec_difference_two[j];
                            best_node.second = dis;
                        }
                    }
                }

            }
            index = 1 - index;
            ++num;
            select_nodes_two.insert(best_node.first);
            cur_nodes[1][best_node.first] = 1;
            hashone[1] += hash_key_temp_one[best_node.first];
            hashtwo[1] += hash_key_temp_two[best_node.first];
            hashthree[1] += hash_key_temp_three[best_node.first];
            Distance max = -1, min = DISTANCE_MAX;
            for (int l = 0; l < nb_nodes; ++l) {
                best_solu_dis_sum[1][l] += ins.dis_nodes(l, best_node.first);
                if (cur_nodes[1][l]) {
                    if (max < best_solu_dis_sum[1][l])max = best_solu_dis_sum[1][l];
                    if (min > best_solu_dis_sum[1][l])min = best_solu_dis_sum[1][l];
                }
            }
            best_max_select_node[1] = max;
            best_min_select_node[1] = min;
            local_best_obj[1] = max - min;
        }
        best_hashfun_one[1] = hashone[0] % size_of_tabu_list;
        best_hashfun_two[1] = hashtwo[0] % size_of_tabu_list;
        best_hashfun_three[1] = hashthree[0] % size_of_tabu_list;
        local_best[1] = cur_nodes[1];
    }

    void CrossOver::tabusearch() {
        //TODO：可以多少步之内改进不了历史最优解停止，目前设置为一定步数
        for (int i = 0; i < population; ++i) {
            bool tabu_flag = false;
            int iter = 0;
            int step_length = 0;         //保存多少步之内改进不了历史最优解
            //while (iter < tabu_iter_max) {
            while(true){
                int _hashfun_one = best_hashfun_one[i];
                int _hashfun_two = best_hashfun_two[i];
                int _hashfun_three = best_hashfun_three[i];
                node_value[i] = local_best[i];               //强化搜索策略：用历史最优解更新当前解
                cur_obj[i] = local_best_obj[i];
                //若历史最优解!=当前解,则要更新当前解的部分数据结构node_dis_sum，max_select_node， min_select_node，select_nodes和no_select_nodes
                max_select_node[i] = best_max_select_node[i];                   //xxf:解决大bug1：之前替换当前解时，未更新相关数据结构node_dis_sum，select_nodes和no_select_nodes
                min_select_node[i] = best_min_select_node[i];
                no_select_nodes[i].clear();            //更新辅助结构select_nodes和no_select_nodes
                select_nodes[i].clear();
                Distance temp_sum = (best_max_select_node[i] + best_min_select_node[i]) / 2;
                for (int j = 0; j < nb_nodes; ++j) {
                    node_dis_sum[i][j] = best_solu_dis_sum[i][j];
                    Distance temp = fabs(best_solu_dis_sum[i][j] - temp_sum);
                    if (node_value[i][j])select_nodes[i].push_back(make_pair(j, temp));
                    else no_select_nodes[i].push_back(make_pair(j, temp));
                }
                sort(no_select_nodes[i].begin(), no_select_nodes[i].end(), compareByAscend);    //未选中的升序排列

                int count = 0;

                while (count <= single_iter_length) {
                    pair<int, int> swap_pair(-1, -1);              //邻域结构中找最好的非禁忌解:保存非禁忌的最好的交换对;第一个I1-->I0，第二个I0-->I1
                    pair<Distance, Distance> new_obj(DISTANCE_MAX, 0);         //保存对应的目标函数的最大距离和最小距离
                    find_best_move(i, swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three);   //_hashfun_one已经是更新解的哈希函数值
                    //TODO:判断是否所有邻域解都在禁忌中
                    if (swap_pair.first == -1) {       //判断是否邻域解都在禁忌中
                        mylog << "邻域解都在禁忌中；当前iter：" << iter <<= logsw_local;
                        tabu_flag = true;
                        break;
                    }
                    if (update_solu(i, swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three, step_length))count = 0;   //更新当前解、历史最优解、count
                    else count++;
                    iter++;
                    if (step_length == tabu_iter_max) {
                        tabu_flag = true;
                        step_length = 0;
                       break;
                    }
                }
                if (tabu_flag)break;
            }

        }
    }

    void CrossOver::update_elite() {
        int min_index = -1;
        if (local_best_obj[0] < local_best_obj[1]) {
            min_index = 0;
        }
        else min_index = 1;
        if (local_best_obj[min_index] < elite_obj[0]) {
            elite[0] = local_best[min_index];
            elite_obj[0] = local_best_obj[min_index];
        }
    }

    void CrossOver::update_best() {
        if (elite_obj[0] < cross_best_obj) {
            cross_best = elite[0];
            cross_best_obj = elite_obj[0];
        }
    }

    void CrossOver::random_elite_one() {
        List<int> node_temp(nb_nodes, 0);
        List<int> select_nodes;
        select_nodes.reserve(nb_sub_nodes);
        int num_select = 0;
        int i = 0;
        time_t t = time(0);
        srand(t);
        while (num_select < nb_sub_nodes) {
            while (true) {
                i = rand() % nb_nodes;
                if (node_temp[i] == 0)break;
            }
            node_temp[i] = 1;
            select_nodes.push_back(i);
            num_select++;
        }
        elite[0] = node_temp;
        Distance min = DISTANCE_MAX, max = 0;
        for (int i = 0; i < nb_sub_nodes; ++i) {
            Distance cur = 0;
            for (int j = 0; j < nb_sub_nodes; ++j) {
                cur += ins.dis_nodes(select_nodes[i], select_nodes[j]);
            }
            if (cur > max)max = cur;
            if (cur < min)min = cur;
        }
        elite_obj[0] = max - min;
    }

    void CrossOver::find_best_move(const int number, pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three) {   //xxf:done,right--12.10
        int new_hashone = _hash_one;     //保存当前解的哈希函数值
        int new_hashtwo = _hash_two;
        int new_hashthree = _hash_three;
        for (int i = 0; i < nb_sub_nodes; ++i)
        {
            int one_toZero_node = select_nodes[number][i].first;
            for (int j = 0; j < size_neighbor_struc; ++j) {
                int zero_toOne_node = no_select_nodes[number][j].first;
                int _new_hash_one = new_hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //计算交换后新解的哈希函数值
                int _new_hash_two = new_hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
                int _new_hash_three = new_hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
                _new_hash_one = (_new_hash_one + size_of_tabu_list) % size_of_tabu_list;                //xxf：解决bug2：防止出现负数和>size_of_tabu_list的数
                _new_hash_two = (_new_hash_two + size_of_tabu_list) % size_of_tabu_list;
                _new_hash_three = (_new_hash_three + size_of_tabu_list) % size_of_tabu_list;
                if (tabu_list_three[_new_hash_three]) {       //xxf：两个或三个哈希函数为了减小冲突碰撞；一个哈希函数容易不同解映射到同一个value值，且很容易出现邻域解都被禁忌情况，因为是从强化搜索策略开始
                    if (tabu_list_two[_new_hash_two])
                        if (tabu_list_one[_new_hash_one]) {
                            continue;
                        }
                }
                Distance temp_min = node_dis_sum[number][zero_toOne_node] - ins.dis_nodes(zero_toOne_node, one_toZero_node);      //记录当前一次动作的最大值和最小值;初始化为交换之后新的选中节点的距离之和
                Distance temp_max = temp_min;
                for (int k = 0; k < nb_sub_nodes; ++k) {
                    if (k == i)continue;
                    int node = select_nodes[number][k].first;
                    Distance update_dis = node_dis_sum[number][node] - ins.dis_nodes(node, one_toZero_node) + ins.dis_nodes(node, zero_toOne_node);
                    if (update_dis > temp_max)temp_max = update_dis;
                    else if (temp_min > update_dis)temp_min = update_dis;
                    else;
                }
                if ((_new_obj.first - _new_obj.second) > (temp_max - temp_min)) {   //更新邻域动作
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
        bool flag = false;         //表示是否更新历史最优解
        node_value[number][_pair.first] = 0;                       //更新当前解
        node_value[number][_pair.second] = 1;
        cur_obj[number] = _new_obj.first - _new_obj.second;
        max_select_node[number] = _new_obj.first;
        min_select_node[number] = _new_obj.second;
        if (local_best_obj[number] > cur_obj[number])       //更新历史最优解，更新历史最优解的相关结构
        {
            local_best[number] = node_value[number];  //如果能改进历史最优解，则更新历史最优解,返回true
            local_best_obj[number] = cur_obj[number];
            best_hashfun_one[number] = _hash_one;
            best_hashfun_two[number] = _hash_two;
            best_hashfun_three[number] = _hash_three;
            flag = true;
            step = 0;
        }
        step++;
        tabu_list_one[_hash_one] = 1;         //更新三个禁忌列表
        tabu_list_two[_hash_two] = 1;
        tabu_list_three[_hash_three] = 1;
        no_select_nodes[number].clear();           //更新辅助结构select_nodes和no_select_nodes
        select_nodes[number].clear();
        Distance temp_sum = (max_select_node[number] + min_select_node[number]) / 2.0;
        if (flag) {           //若更新了历史最优解，则保存历史最优解的相关数据结构best_solu_dis_sum
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
        sort(no_select_nodes[number].begin(), no_select_nodes[number].end(), compareByAscend);    //未选中的升序排列
        return flag;
    }
}