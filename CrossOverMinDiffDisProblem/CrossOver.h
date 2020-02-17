#pragma once
#ifndef MIN_DIFF_DP_LOCALSEARCH_H
#define MIN_DIFF_DP_LOACLSEARCH_H
#include "InputOutput.h"

using namespace std;

namespace min_diff_dp {

    class CrossOver {
    public:
        CrossOver(const UMatrix &_matrix, int _tabu_iter, int _tabu_length, int _iter_cycle, int popu, double param, const Solution &_init_sol_one, const Solution &_init_sol_two,
            const Solution &_init_elite_one, const Solution &_init_elite_two, const int _size_of_tabu, const double _param1, const double _param2, const double _param3); 
        Solution solve();
    private:
        void init();
        void cross();
        void tabusearch();
        void update_elite();   //更新elite1解
        void update_best();      //跟新历史最优解
        void random_elite_one();
    private:
        void find_best_move(const int number, pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three);
        bool update_solu(const int number, const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, int &step);

    private:
        const UMatrix &ins;      //算例矩阵

        List<int> cross_best;
        Distance cross_best_obj;

        List<List<int>> elite;
        List<Distance> elite_obj;

        int iter_cycle;
        int cycle;

        List<int> tabu_list_one;                   //三个禁忌列表
        List<int> tabu_list_two;
        List<int> tabu_list_three;
        List<int> hash_key_temp_one;          //保存所有节点的hash中间键值(int)(floor(pow(i, hashFun_one_param)))
        List<int> hash_key_temp_two;
        List<int> hash_key_temp_three;

        int nb_nodes;                   //图中节点数目
        int nb_sub_nodes;               //所选的节点数目

        int population;               //种群个数
        int single_iter_length;                      //禁忌搜索长度
        int tabu_iter_max;                     //禁忌搜索的最大迭代次数
        int generation;                      //交叉次数
        int size_neighbor_struc;          //参与邻域动作的I0元素个数

        long long max_time;             //迭代运行的最长时间=算例的节点数目
        double rate_of_sele_nodes;        //从I0中选择邻域动作的大小比例
        int size_of_tabu_list;        //L大小
        double hashFun_one_param;     //三个哈希函数的参数
        double hashFun_two_param;
        double hashFun_three_param;

        List<List<int>> node_value;    //当前解        ？（是用List保存，还是Set保存；更快）
        List<Distance> cur_obj;        //当前目标
        List<List<int>> local_best;     //局部最优解
        List<Distance> local_best_obj;        //局部最优目标函数

        List<List<Distance>> node_dis_sum;     //每个节点与选中集合中节点的距离和--邻域结构
        List<Distance> max_select_node;                  //记录node_dis_sum中距离最大值的选中节点
        List<Distance> min_select_node;                  //记录node_dis_sum中距离最小值的选中节点
                                     //保存历史最优解的相关数据；强化搜索策略（S<-历史最优解）中方便使用
        List<List<Distance>> best_solu_dis_sum;
        List<Distance> best_max_select_node;
        List<Distance> best_min_select_node;

        List<List<Pair<int, Distance>>> no_select_nodes;         //未选中的节点排序--辅助结构,Dis-(max+min)/2
        List<List<Pair<int, Distance>>> select_nodes;             //选中的节点排序--辅助结构

        List<int> best_hashfun_one;      //中间值:历史最优解的三个哈希函数值
        List<int> best_hashfun_two;
        List<int> best_hashfun_three;

    };

}

#endif