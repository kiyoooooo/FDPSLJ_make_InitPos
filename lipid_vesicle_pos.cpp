#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <random>
#include <string.h>//文字列の代入に使う
#include <bits/stdc++.h>//piの利用で必要(M_PI)


int main() {

  //自分で決めるパラメータ
  const int    box_x_y = 1 ;//ボックスサイズのxが1としたときのyの比
  const int    box_x_z = 1 ;//ボックスサイズのxが1としたときのzの比
  const int    num_inner_water = 1200000  ;//脂質の内側に存在する水粒子の数
  //1100 11000 110000 1200000
  const int    num_outer_water = 25000000 ;//脂質の外側に存在する水粒子の数
  //24000 240000 2400000 25000000
  const int    num_inner_lipid_part  = 440  ;//ベシクルの内側を構成する脂質の2次元半円に含まれる数 約　*0.32　で下の1/10のサイズ
  //15 45 140 440
  const int    num_outer_lipid_part  = 530  ;//ベシクルの外側を構成する脂質の2次元半円に含まれる数
  //17 55 170 530
  const double bond_length     = 0.5 ;//脂質分子間の長さ
  const double lipid_to_lipid_length = 0.5 ;//内脂質分子と外脂質分子の距離
  const double vesicle_radius  = 13.0* 2.154 * 2.154 ;// /2.154 ; // * 2.154 * 2.154 ;//ベシクルの半径(inner_lipidとouter_lipidの中間と原点を測る) 2.154 ~= powf(10,1/3) 粒子数10倍にしたければ10倍の(^1/3)をかけたもの．
  const double vesicle_core_pos_x = 0;//ベシクル中心x座標
  const double vesicle_core_pos_y = 0;//ベシクル中心y座標
  const double vesicle_core_pos_z = 0;//ベシクル中心z座標



  
  //確定するパラメータ
  const int    num_inner_lipid = 2 * (num_inner_lipid_part - 1) * (num_inner_lipid_part - 2) + 2;//ベシクルの内側を構成する脂質の数
  const int    num_outer_lipid = 2 * (num_outer_lipid_part - 1) * (num_outer_lipid_part - 2) + 2;//ベシクルの外側を構成する脂質の数
  const int    num_water = num_inner_water + num_outer_water;//水粒子の数
  const int    num_lipid = (num_inner_lipid + num_outer_lipid) * 4;//脂質を構成する粒子の数
  const int    num = num_inner_water + num_outer_water + (num_inner_lipid + num_outer_lipid) * 4;//すべての粒子数
  const double rho = 3.0;//密度
  const double box_size_x = powf((double)num / rho, 1.0 / 3.0) / powf(box_x_y * box_x_z, 1.0 / 3.0);//ボックスサイズxyz三辺が1：1：1の場合．今ｘが1なので，box_size_xになる．
  const double box_size_y = box_size_x * box_x_y;
  const double box_size_z = box_size_x * box_x_z;
  const double vesicle_thickness = bond_length * 6 + lipid_to_lipid_length;//ベシクルを構成する脂質二重層の厚み
  const double vesicle_inner_radius2 = powf(vesicle_radius - lipid_to_lipid_length / 2 - bond_length * 3.0, 2.0);//ベシクルの半径(ベシクルの中心とinner_lipidまでの距離の2乗)
  const double vesicle_outer_radius2 = powf(vesicle_radius + lipid_to_lipid_length / 2 + bond_length * 3.0, 2.0);//ベシクルの半径(ベシクルの中心とouter_lipidまでの距離の2乗)
  const double angle_to_place_inner_lipid = M_PI / (num_inner_lipid_part - 1);//ベシクルをの内側を構成する脂質の角度間隔
  const double angle_to_place_outer_lipid = M_PI / (num_outer_lipid_part - 1);//ベシクルをの外側を構成する脂質の角度間隔


  std::string filename0 = "../output/detail.txt";
  std::ofstream writing_file0;
  writing_file0.open(filename0, std::ios::out);




  writing_file0 << "自分で決めるパラメータ" << std::endl;
  writing_file0 << "box_x_y = " << box_x_y << " (ボックスサイズのxが1としたときのyの比) " << std::endl;
  writing_file0 << "box_x_z = " << box_x_z << " (ボックスサイズのxが1としたときのzの比) " << std::endl;
  writing_file0 << "num_inner_water = " << num_inner_water << " (脂質の内側に存在する水粒子の数) " << std::endl;
  writing_file0 << "num_outer_water = " << num_outer_water << " (脂質の外側に存在する水粒子の数) " << std::endl;
  writing_file0 << "num_inner_lipid_part = " << num_inner_lipid_part << " (ベシクルの内側を構成する脂質の2次元半円に含まれる数) " << std::endl;
  writing_file0 << "num_outer_lipid_part = " << num_outer_lipid_part << " (ベシクルの外側を構成する脂質の2次元半円に含まれる数) " << std::endl;
  writing_file0 << "bond_length = " << bond_length << " (脂質分子間の長さ) " << std::endl;
  writing_file0 << "lipid_to_lipid_length = " << lipid_to_lipid_length << " (内脂質分子と外脂質分子の距離) " << std::endl;
  writing_file0 << "vesicle_radius = " << vesicle_radius << " (ベシクルの半径(inner_lipidとouter_lipidの中間と原点を測る) 2.154 ~= powf(10,1/3)) " << std::endl;
  writing_file0 << "vesicle_core_pos_x = " << vesicle_core_pos_x << " (ベシクル中心x座標) " << std::endl;
  writing_file0 << "vesicle_core_pos_y = " << vesicle_core_pos_y << " (ベシクル中心y座標) " << std::endl;
  writing_file0 << "vesicle_core_pos_z = " << vesicle_core_pos_z << " (ベシクル中心z座標) " << std::endl;

  writing_file0 << std::endl;

  writing_file0 << "確定するパラメータ" << std::endl;
  writing_file0 << "num_inner_lipid = " << num_inner_lipid << " (ベシクルの内側を構成する脂質の数) " << std::endl;
  writing_file0 << "num_outer_lipid = " << num_outer_lipid << " (ベシクルの外側を構成する脂質の数) " << std::endl;
  writing_file0 << "num_water = " << num_water << " (水粒子の数) " << std::endl;
  writing_file0 << "num_lipid = " << num_lipid << " (脂質を構成する粒子の数) " << std::endl;
  writing_file0 << "num = " << num << " (すべての粒子数) " << std::endl;
  writing_file0 << "rho = " << rho << " (密度) " << std::endl;
  writing_file0 << "box_size_x = " << box_size_x << " (ボックスサイズxyz三辺が1：1：1の場合, box_size_x = box_size_y = box_size_z. ) " << std::endl;
  writing_file0 << "box_size_y = " << box_size_y << std::endl;
  writing_file0 << "box_size_z = " << box_size_z << std::endl;
  writing_file0 << "vesicle_thickness = " << vesicle_thickness << " (ベシクルを構成する脂質二重層の厚み) " << std::endl;
  writing_file0 << "vesicle_inner_radius2 = " << vesicle_inner_radius2 << " (ベシクルの半径(ベシクルの中心とinner_lipidまでの距離の2乗)) " << std::endl;
  writing_file0 << "vesicle_outer_radius2 = " << vesicle_outer_radius2 << " (ベシクルの半径(ベシクルの中心とouter_lipidまでの距離の2乗)) " << std::endl;
  writing_file0 << "angle_to_place_inner_lipid = " << angle_to_place_inner_lipid << " (ベシクルをの内側を構成する脂質の角度間隔) " << std::endl;
  writing_file0 << "angle_to_place_outer_lipid = " << angle_to_place_outer_lipid << " (ベシクルをの外側を構成する脂質の角度間隔) " << std::endl;

  writing_file0 << std::endl;
  
  writing_file0 << "water : lipid = " << (double)num_water / (double)num * 100.0 << " : " << (double)num_lipid / (double)num * 100.0 << std::endl;

  
  if(box_size_x < std::sqrt(vesicle_outer_radius2) * 2) {
    std::cout << "NG   box_size_x = " << box_size_x << "   sqrt(vesicle_outer_radius2) * 2 = " << sqrt(vesicle_outer_radius2) * 2 << std::endl;
    return 0;
  }

  else {
    std::cout << "OK   box_size_x = " << box_size_x << "   sqrt(vesicle_outer_radius2) * 2 = " << sqrt(vesicle_outer_radius2) * 2 << std::endl;
    writing_file0 << "OK   box_size_x = " << box_size_x << "   >   sqrt(vesicle_outer_radius2) * 2 = " << sqrt(vesicle_outer_radius2) * 2 << std::endl;
  }
  
  //open file
  //pos_file
  FILE *fpo0;
  fpo0 = fopen("../output/initial_pos_lipid.cdv", "a");
  if(fpo0 == NULL) {
    printf("ERROR_initial_pos_lipid.cdv\n");
    return -1;
  }


  //vel_file
  FILE *fpo1;
  fpo1 = fopen("../output/initial_vel_lipid.cdv", "a");
  if(fpo1 == NULL) {
    printf("ERROR_initial_vel_lipid.cdv\n");
    return -1;
  }


  //bond_file
   FILE *fpo2;
  fpo2 = fopen("../output/bond_info.cdv", "a");
  if(fpo2 == NULL) {
    printf("ERROR_bond_info.cdv\n");
    return -1;
  }


  //angle_file
    FILE *fpo3;
  fpo3 = fopen("../output/angle_info.cdv", "a");
  if(fpo3 == NULL) {
    printf("ERROR_angle_info.cdv\n");
    return -1;
  }



  
  
  //fprintf(fpo0, "hello\n");
  //fclose(fpo0);
  fprintf(fpo0, "'box_sx=%lf box_sy=%lf box_sz=%lf box_ex=%lf box_ey=%lf box_ez=%lf box_wt=.05\n", - box_size_x / 2 , - box_size_y / 2, - box_size_z / 2, box_size_x / 2 , box_size_y / 2, box_size_z / 2);
  fprintf(fpo0, "'r1=0.4 r2=0.4 r3=0.\n");
  fprintf(fpo0, "'c1=(0.5,0.5,1) c2=(1,0.5,0.5) c3=(0,1,1)\n");



  
  int particle_id       = 1;//粒子番号
  int type_lipid_a      = 1;//粒子種，脂質親水基
  int type_lipid_b      = 2;//粒子種，脂質疎水基
  int type_inner_water  = 3;//粒子種，内側の水
  int type_outer_water  = 3;//粒子種，外側の水
  int type_bond_a       = 1;//lipidのbondの種類（現状1種）
  int type_angle_a      = 1;//lipidのangleの種類（現状1種）
  double temp_angle_xy;// (= angle_to_place_inner_lipid)
  double temp_angle_xz;// (= angle_to_place_inner_lipid)
  double lipid_delete_rate = 0.4;
  double size_ratio = 1.0/1.7;
  double power_ratio = 10.0;
  std::random_device seed_gen;
  std::mt19937_64 engine(seed_gen());
  // [0.0, 1.0) の一様分布実数生成器
  std::uniform_real_distribution<double> get_rand_uni_real(0.0, 1.0);



  //lipidの配置
  //ベシクル内側を構成するlipid

  for(temp_angle_xy = 0; temp_angle_xy < M_PI; temp_angle_xy += angle_to_place_inner_lipid) {
    for(temp_angle_xz = 0; temp_angle_xz < 2 * M_PI; temp_angle_xz += angle_to_place_inner_lipid) {


//std::cout<<(1/2)*(1-std::cos(temp_angle_xz))<<std::endl;
//std::cout<<(1.0/2.0)*(1-std::cos(temp_angle_xz))<<std::endl;// 1.0が必要だった．そうでなければ，整数型にキャストされてゼロになってしまう．
//std::cout<<temp_angle_xz<<std::endl;
//std::cout<<(1/2)*(1-std::cos(2*M_PI*2*temp_angle_xz/2*M_PI))* get_rand_uni_real(engine)<<std::endl;
//std::cout<< get_rand_uni_real(engine)<<std::endl;

//      if(temp_angle_xz *lipid_delete_rate*get_rand_uni_real(engine) > 2 * M_PI) {//密度の高い部分を改善するため
//      if(temp_angle_xz <=  angle_to_place_outer_lipid*10) {//密度の高い部分を改善するため
//      if(std::pow((1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz)),8.0)* get_rand_uni_real(engine) <=  lipid_delete_rate && (1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz)) != 1.0) {//密度の高い部分を改善するため
//      if(-(std::pow((1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz-M_PI)),2.0)-1)* get_rand_uni_real(engine) <=  lipid_delete_rate && (1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz)) != 1.0) {//密度の高い部分を改善するため
        if(std::pow(size_ratio*(1.0-std::cos(2.0*temp_angle_xz)),power_ratio)* get_rand_uni_real(engine) <=  lipid_delete_rate && (1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz)) != 1.0) {//密度の高い部分を改善するため
      //疎水基3_write_pos
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id, type_lipid_b, \
	      (vesicle_core_pos_x + vesicle_radius - lipid_to_lipid_length / 2) * std::cos(temp_angle_xz) * std::cos(temp_angle_xy), \
	      (vesicle_core_pos_y + vesicle_radius - lipid_to_lipid_length / 2) * std::sin(temp_angle_xy) * std::cos(temp_angle_xz), \
	      (vesicle_core_pos_z + vesicle_radius - lipid_to_lipid_length / 2) * std::sin(temp_angle_xz));
      //疎水基3_write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id, 0.0, 0.0, 0.0);



      //疎水基2_write_pos
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id + 1, type_lipid_b, \
	      (vesicle_core_pos_x + vesicle_radius - lipid_to_lipid_length / 2 - bond_length) * std::cos(temp_angle_xz) * std::cos(temp_angle_xy), \
	      (vesicle_core_pos_y + vesicle_radius - lipid_to_lipid_length / 2 - bond_length) * std::sin(temp_angle_xy) * std::cos(temp_angle_xz), \
	      (vesicle_core_pos_z + vesicle_radius - lipid_to_lipid_length / 2 - bond_length) * std::sin(temp_angle_xz));
      //疎水基2_write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id + 1, 0.0, 0.0, 0.0);
      //疎水基3＋疎水基2_write_bond
      fprintf(fpo2, "%d   %d   %d \n", particle_id, particle_id + 1, type_bond_a);

      

      //疎水基1_write_pos 
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id + 2, type_lipid_b, \
	      (vesicle_core_pos_x + vesicle_radius - lipid_to_lipid_length / 2 - bond_length * 2) * std::cos(temp_angle_xz) * std::cos(temp_angle_xy), \
	      (vesicle_core_pos_y + vesicle_radius - lipid_to_lipid_length / 2 - bond_length * 2) * std::sin(temp_angle_xy) * std::cos(temp_angle_xz), \
	      (vesicle_core_pos_z + vesicle_radius - lipid_to_lipid_length / 2 - bond_length * 2) * std::sin(temp_angle_xz));
      //疎水基1_write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id + 2, 0.0, 0.0, 0.0);
      //疎水基2＋疎水基1_write_bond
      fprintf(fpo2, "%d   %d   %d \n", particle_id + 1, particle_id + 2, type_bond_a);
      //疎水基3＋疎水基2＋疎水基1_write_angle
      fprintf(fpo3, "%d   %d   %d   %d \n", particle_id, particle_id + 1, particle_id + 2, type_angle_a);

      
      //親水基_write_pos
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id + 3, type_lipid_a, \
	      (vesicle_core_pos_x + vesicle_radius - lipid_to_lipid_length / 2 - bond_length * 3) * std::cos(temp_angle_xz) * std::cos(temp_angle_xy), \
	      (vesicle_core_pos_y + vesicle_radius - lipid_to_lipid_length / 2 - bond_length * 3) * std::sin(temp_angle_xy) * std::cos(temp_angle_xz), \
	      (vesicle_core_pos_z + vesicle_radius - lipid_to_lipid_length / 2 - bond_length * 3) * std::sin(temp_angle_xz));
      //親水基_write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id + 3, 0.0, 0.0, 0.0);
      //疎水基1＋親水基_write_bond
      fprintf(fpo2, "%d   %d   %d \n", particle_id + 2, particle_id + 3, type_bond_a);
      //疎水基2＋疎水基1＋親水基_write_angle
      fprintf(fpo3, "%d   %d   %d   %d \n", particle_id + 1, particle_id + 2, particle_id + 3, type_angle_a);

      particle_id += 4;
      }
    }
  }

  //ベシクルの外側を構成するlipid                                                                                                                                                                                                     

  for(temp_angle_xy = 0; temp_angle_xy < M_PI; temp_angle_xy += angle_to_place_outer_lipid) {
    for(temp_angle_xz = 0; temp_angle_xz < 2 * M_PI; temp_angle_xz += angle_to_place_outer_lipid) {


//      if(temp_angle_xz *lipid_delete_rate*get_rand_uni_real(engine) > 2 * M_PI) {//密度の高い部分を改善するため
//      if(temp_angle_xz <=  angle_to_place_outer_lipid*10) {//密度の高い部分を改善するため
//      if((1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz))* get_rand_uni_real(engine) <=  lipid_delete_rate) {//密度の高い部分を改善するため
//      if(-(std::pow((1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz-M_PI)),2.0)-1)* get_rand_uni_real(engine) <=  lipid_delete_rate && (1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz)) != 1.0) {//密度の高い部分を改善するため
        if(std::pow(size_ratio*(1.0-std::cos(2.0*temp_angle_xz)),power_ratio)* get_rand_uni_real(engine) <=  lipid_delete_rate && (1.0/2.0)*(1.0-std::cos(2.0*temp_angle_xz)) != 1.0) {//密度の高い部分を改善するため
      //疎水基3_write_pos                                                                                                                                                                                    
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id, type_lipid_b, \
              (vesicle_core_pos_x + vesicle_radius + lipid_to_lipid_length / 2) * std::cos(temp_angle_xz) * std::cos(temp_angle_xy), \
              (vesicle_core_pos_y + vesicle_radius + lipid_to_lipid_length / 2) * std::sin(temp_angle_xy) * std::cos(temp_angle_xz), \
              (vesicle_core_pos_z + vesicle_radius + lipid_to_lipid_length / 2) * std::sin(temp_angle_xz));
      //疎水基3_write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id, 0.0, 0.0, 0.0);

      //疎水基2_write_pos                                                                                                                                                                                   
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id + 1, type_lipid_b, \
              (vesicle_core_pos_x + vesicle_radius + lipid_to_lipid_length / 2 + bond_length) * std::cos(temp_angle_xz) * std::cos(temp_angle_xy), \
              (vesicle_core_pos_y + vesicle_radius + lipid_to_lipid_length / 2 + bond_length) * std::sin(temp_angle_xy) * std::cos(temp_angle_xz), \
              (vesicle_core_pos_z + vesicle_radius + lipid_to_lipid_length / 2 + bond_length) * std::sin(temp_angle_xz));
      //疎水基2_write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id + 1, 0.0, 0.0, 0.0);
      //疎水基3＋疎水基2_write_bond
      fprintf(fpo2, "%d   %d   %d \n", particle_id, particle_id + 1, type_bond_a);


      //疎水基1_write_pos                                                                                                                                                                                             
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id + 2, type_lipid_b, \
              (vesicle_core_pos_x + vesicle_radius + lipid_to_lipid_length / 2 + bond_length * 2) * std::cos(temp_angle_xz) * std::cos(temp_angle_xy), \
              (vesicle_core_pos_y + vesicle_radius + lipid_to_lipid_length / 2 + bond_length * 2) * std::sin(temp_angle_xy) * std::cos(temp_angle_xz), \
              (vesicle_core_pos_z + vesicle_radius + lipid_to_lipid_length / 2 + bond_length * 2) * std::sin(temp_angle_xz));
      //疎水基1_write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id + 2, 0.0, 0.0, 0.0);
      //疎水基2＋疎水基1_write_bond
      fprintf(fpo2, "%d   %d   %d \n", particle_id + 1, particle_id + 2, type_bond_a);
      //疎水基3＋疎水基2＋疎水基1_write_angle
      fprintf(fpo3, "%d   %d   %d   %d \n", particle_id, particle_id + 1, particle_id + 2, type_angle_a);



      //親水基_write_pos                                                                                                                                                                                              
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id + 3, type_lipid_a, \
              (vesicle_core_pos_x + vesicle_radius + lipid_to_lipid_length / 2 + bond_length * 3) * std::cos(temp_angle_xz) * std::cos(temp_angle_xy), \
              (vesicle_core_pos_y + vesicle_radius + lipid_to_lipid_length / 2 + bond_length * 3) * std::sin(temp_angle_xy) * std::cos(temp_angle_xz), \
              (vesicle_core_pos_z + vesicle_radius + lipid_to_lipid_length / 2 + bond_length * 3) * std::sin(temp_angle_xz));
      //親水基_write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id + 3, 0.0, 0.0, 0.0);
      //疎水基1＋親水基_write_bond
      fprintf(fpo2, "%d   %d   %d \n", particle_id + 2, particle_id + 3, type_bond_a);
      //疎水基2＋疎水基1＋親水基_write_angle
      fprintf(fpo3, "%d   %d   %d   %d \n", particle_id + 1, particle_id + 2, particle_id + 3, type_angle_a);

      particle_id += 4;
      }
    }
  }

  //水の配置
  
  int temp_num_inner_water = 0;
  int temp_num_outer_water = 0;
  double random_x;
  double random_y;
  double random_z;
  //上で定義したため，ここではコメント化
//  std::random_device seed_gen;
//  std::mt19937_64 engine(seed_gen());
  // [0.0, 1.0) の一様分布実数生成器
//  std::uniform_real_distribution<double> get_rand_uni_real(0.0, 1.0);

    

  while(temp_num_inner_water < num_inner_water || temp_num_outer_water < num_outer_water) {



    random_x = (get_rand_uni_real(engine) - 0.50) * box_size_x;//本来は (get_rand_uni_real(engine) - 0.50) * 2 * box_size_x / 2
    random_y = (get_rand_uni_real(engine) - 0.50) * box_size_y;
    random_z = (get_rand_uni_real(engine) - 0.50) * box_size_z;//ボックスの中心が原点(0,0)であることに注意．負の座標が許されないなら座標をbox_size/2ずらす必要がある．


    
    //ベシクルの内側に存在する水粒子
    if(powf(random_x - vesicle_core_pos_x, 2.0) + powf(random_y - vesicle_core_pos_y, 2.0) + powf(random_z - vesicle_core_pos_z, 2.0) < vesicle_inner_radius2 && temp_num_inner_water < num_inner_water) {
      //write_pos
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id, type_inner_water, \
	      random_x, \
	      random_y, \
	      random_z);

      //write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id, 0.0, 0.0, 0.0);

      particle_id += 1;
      temp_num_inner_water++;
    }
    //ベシクルの外側に存在する水粒子
    else if(powf(random_x - vesicle_core_pos_x, 2.0) + powf(random_y - vesicle_core_pos_y, 2.0) + powf(random_z - vesicle_core_pos_z, 2.0) > vesicle_outer_radius2 && temp_num_outer_water < num_outer_water) {
      //_write_pos
      fprintf(fpo0, "%d %d   %lf   %lf   %lf \n", particle_id, type_outer_water, \
              random_x, \
              random_y, \
              random_z);

      //write_vel
      fprintf(fpo1, "%d   %lf   %lf   %lf \n", particle_id, 0.0, 0.0, 0.0);

      particle_id += 1;
      temp_num_outer_water++;
    }

    //    std::cout<<"temp_num_inner_water"<<temp_num_inner_water<<"temp_num_outer_water"<<temp_num_outer_water<<std::endl;
  }










  












#ifdef UNI_SET
  double unit_size = cell_size/(double)nunit;
  double ush = unit_size * 0.5;
  PS::F64vec unit[4];
  unit[0].x = 0.0; unit[1].x = ush; unit[2].x = 0.0; unit[3].x = ush;
  unit[0].y = 0.0; unit[1].y = ush; unit[2].y = ush; unit[3].y = 0.0;
  unit[0].z = 0.0; unit[1].z = 0.0; unit[2].z = ush; unit[3].z = ush;

  int ip=0;
  for(int i=0; i<nunit; i++){
    for(int j=0; j<nunit; j++){
      for(int k=0; k<nunit; k++){
        for(int l=0; l<4; l++){
          pos[ip].x = i*unit_size + 0.1*unit[l].x + 0.1*ush;
          pos[ip].y = j*unit_size + 0.1*unit[l].y + 0.1*ush;
          pos[ip].z = k*unit_size + 0.1*unit[l].z + 0.1*ush;
          ip++;
        }
      }
    }
  }
#endif


#ifdef RANDOM_SET
  PS::F64 random_x = 0;
  PS::F64 random_y = 0;
  PS::F64 random_z = 0;


  int ip=0;
  for(int i=0; i*12<n_tot;i++){
    random_x = (mt.genrand_res53() - 0.50) * cell_size;
    random_y = (mt.genrand_res53() - 0.50) * cell_size;
    random_z = (mt.genrand_res53() - 0.50) * cell_size;
    for(int l=0; l<12; l++){
      pos[ip].x = random_x+l*0.80;
      pos[ip].y = random_y;
      pos[ip].z = random_z;
      ip++;
    }
  }
#endif


return 0; 
}
