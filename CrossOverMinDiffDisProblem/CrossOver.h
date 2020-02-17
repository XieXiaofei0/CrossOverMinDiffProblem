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
        void update_elite();   //����elite1��
        void update_best();      //������ʷ���Ž�
        void random_elite_one();
    private:
        void find_best_move(const int number, pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three);
        bool update_solu(const int number, const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, int &step);

    private:
        const UMatrix &ins;      //��������

        List<int> cross_best;
        Distance cross_best_obj;

        List<List<int>> elite;
        List<Distance> elite_obj;

        int iter_cycle;
        int cycle;

        List<int> tabu_list_one;                   //���������б�
        List<int> tabu_list_two;
        List<int> tabu_list_three;
        List<int> hash_key_temp_one;          //�������нڵ��hash�м��ֵ(int)(floor(pow(i, hashFun_one_param)))
        List<int> hash_key_temp_two;
        List<int> hash_key_temp_three;

        int nb_nodes;                   //ͼ�нڵ���Ŀ
        int nb_sub_nodes;               //��ѡ�Ľڵ���Ŀ

        int population;               //��Ⱥ����
        int single_iter_length;                      //������������
        int tabu_iter_max;                     //��������������������
        int generation;                      //�������
        int size_neighbor_struc;          //������������I0Ԫ�ظ���

        long long max_time;             //�������е��ʱ��=�����Ľڵ���Ŀ
        double rate_of_sele_nodes;        //��I0��ѡ���������Ĵ�С����
        int size_of_tabu_list;        //L��С
        double hashFun_one_param;     //������ϣ�����Ĳ���
        double hashFun_two_param;
        double hashFun_three_param;

        List<List<int>> node_value;    //��ǰ��        ��������List���棬����Set���棻���죩
        List<Distance> cur_obj;        //��ǰĿ��
        List<List<int>> local_best;     //�ֲ����Ž�
        List<Distance> local_best_obj;        //�ֲ�����Ŀ�꺯��

        List<List<Distance>> node_dis_sum;     //ÿ���ڵ���ѡ�м����нڵ�ľ����--����ṹ
        List<Distance> max_select_node;                  //��¼node_dis_sum�о������ֵ��ѡ�нڵ�
        List<Distance> min_select_node;                  //��¼node_dis_sum�о�����Сֵ��ѡ�нڵ�
                                     //������ʷ���Ž��������ݣ�ǿ���������ԣ�S<-��ʷ���Ž⣩�з���ʹ��
        List<List<Distance>> best_solu_dis_sum;
        List<Distance> best_max_select_node;
        List<Distance> best_min_select_node;

        List<List<Pair<int, Distance>>> no_select_nodes;         //δѡ�еĽڵ�����--�����ṹ,Dis-(max+min)/2
        List<List<Pair<int, Distance>>> select_nodes;             //ѡ�еĽڵ�����--�����ṹ

        List<int> best_hashfun_one;      //�м�ֵ:��ʷ���Ž��������ϣ����ֵ
        List<int> best_hashfun_two;
        List<int> best_hashfun_three;

    };

}

#endif