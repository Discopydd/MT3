#pragma once
#include "Vector3.h"
#include "Matrix4x4.h"
struct Line {
	Vector3 origir;//!<始点
	Vector3 diff;//!<終点への差分ベクトル
};

struct Ray {
	Vector3 origin;//!<始点
	Vector3 diff;//!<終点への差分ベクトル 
};

struct Segment {
	Vector3 origin;//!<始点
	Vector3 diff;//!<終点への差分ベクトル
};

struct Sphere {
    Vector3 center;
    float radius;
};

struct Plane
{
	Vector3 normal;//!<法線
	float distance;//!<距離
};
//三角形
struct Triangle {
    Vector3 vertices[3]; //!< 顶点
};
struct AABB {
	Vector3 min;
	Vector3 max;

};

struct Spring {
    Vector3 anchor;      // 固定された端の位置
    float naturalLength; // 自然長
    float stiffness;     // 剛性. バネ定数k
    float dampingCoefficient; // 減衰係数
};

struct Ball {
    Vector3 position;    // ボールの位置
    Vector3 velocity;    // ボールの速度
    Vector3 acceleration;// ボールの加速度
    float mass;          // ボールの質量
    float radius;        // ボールの半径
    unsigned int color;  // ボールの色
};
struct Pendulum {
    Vector3 anchor;             // 固定された端の位置
    float length;               // 紐の長さ
    float angle;                // 現在の角度
    float angularVelocity;      // 角速度
    float angularAcceleration;  // 角加速度
};