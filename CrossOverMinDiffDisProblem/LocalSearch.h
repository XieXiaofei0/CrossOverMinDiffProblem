#pragma once
#ifndef MIN_DIFF_DP_LOCALSEARCH_H
#define MIN_DIFF_DP_LOACLSEARCH_H
#include "InputOutput.h"

using namespace std;

namespace min_diff_dp{

class LocalSearch {

public:
    LocalSearch(const UMatrix &_matrix, double param, const Solution &_init_sol, const int _tabu_step,
        const int _size_of_tabu, const double _param1, const double _param2, const double _param3);
    Solution solve();
    Solution solve_tabu();

private:
    void init();         //初始化数据结构

private:
    void find_best_move(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three);   //局部搜索中在邻域结构中找到最好的交换动作
    bool update_solu(const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, int &step);   //做完邻域动作后更新当前解和历史最优解
    //加入点的禁忌后的查找和更新
    bool find_best_move_tabu(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, const int &iter);
    bool update_solu_tabu(bool tabu_flag, const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, int &step, int tabu_length);
    //test
    //bool stochastic_perturbation(int &hashone, int &hashtwo, int &hashthree, Distance &max, Distance &min);
    //test end

private:
    int hash_function_one();
    int hash_function_two();
    int hash_function_three();
//    int hash_function_temp_one(const vector<int>& temp);
//    int hash_function_temp_two(const vector<int>& temp);
//    int hash_function_temp_three(const vector<int>& temp);

private:
    const UMatrix &ins;      //算例矩阵

    List<int> node_value;    //当前解        ？（是用List保存，还是Set保存；更快）
    Distance cur_obj;        //当前目标
    List<int> local_best;     //局部最优解
    Distance local_best_obj;        //局部最优目标函数

    List<Distance> node_dis_sum;     //每个节点与选中集合中节点的距离和--邻域结构
    Distance max_select_node;                  //记录node_dis_sum中距离最大值的选中节点
    Distance min_select_node;                  //记录node_dis_sum中距离最小值的选中节点
                                 //保存历史最优解的相关数据；强化搜索策略（S<-历史最优解）中方便使用
    List<Distance> best_solu_dis_sum;
    Distance best_max_select_node;
    Distance best_min_select_node;

    List<Pair<int, Distance>> no_select_nodes;         //未选中的节点排序--辅助结构,Dis-(max+min)/2
    List<Pair<int, Distance>> select_nodes;             //选中的节点排序--辅助结构
    List<int> tabu_list_one;                   //三个禁忌列表
    List<int> tabu_list_two;
    List<int> tabu_list_three;
    List<int> hash_key_temp_one;          //保存所有节点的hash中间键值(int)(floor(pow(i, hashFun_one_param)))
    List<int> hash_key_temp_two;
    List<int> hash_key_temp_three;

    //点的禁忌
    List<int> tabu_nodes;                

    int best_hashfun_one;      //中间值:历史最优解的三个哈希函数值
    int best_hashfun_two;
    int best_hashfun_three;   

    int nb_nodes;                   //图中节点数目
    int nb_sub_nodes;               //所选的节点数目

    int iter;                      //总迭代次数
    int size_neighbor_struc;          //参与邻域动作的I0元素个数

    //参数相关
    long long max_time;             //迭代运行的最长时间=算例的节点数目
    double rate_of_sele_nodes;        //从I0中选择邻域动作的大小比例
    int tabu_step;                 //禁忌步长
    int size_of_tabu_list;        //L大小
    double hashFun_one_param;     //三个哈希函数的参数
    double hashFun_two_param;
    double hashFun_three_param;
};

}

#endif