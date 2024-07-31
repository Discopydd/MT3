#include <Novice.h>
#include<imgui.h>
#include"MyMath.h"
#include "Vector3.h"
#include "Matrix4x4.h"


static const int kWindowWidth = 1280; 
static const int kWindowHeight = 720;

const char kWindowTitle[] = "LE2C_29_リ_ヨン";


bool isStarted = false; 

Spring spring{};
Ball ball{};

void Reset() {
    spring.anchor = {0.0f, 0.6f, 0.0f};
    spring.naturalLength = 1.0f;
    spring.stiffness = 100.0f;
    spring.dampingCoefficient = 2.0f; 

    ball.position = {1.2f, 0.6f, 0.0f};
    ball.velocity = {0.0f, 0.0f, 0.0f};
    ball.acceleration = {0.0f, 0.0f, 0.0f};
    ball.mass = 2.0f;
    ball.radius = 0.05f;
    ball.color = 0x0000FFFF;

    isStarted = false;
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

    // ライブラリの初期化
    Novice::Initialize(kWindowTitle, 1280, 720);


    Vector3 cameraTranslate = {0.0f, 2.64f, -9.27f};
    Vector3 cameraRotate = {0.26f, 0.0f, 0.0f};



    Reset(); 


    float deltaTime = 1.0f / 60.0f;

    // キー入力結果を受け取る箱
    char keys[256] = {0};
    char preKeys[256] = {0};

    // ウィンドウの×ボタンが押されるまでループ
    while (Novice::ProcessMessage() == 0) {
        // フレームの开始
        Novice::BeginFrame();

        // キー入力を受け取る
        memcpy(preKeys, keys, 256);
        Novice::GetHitKeyStateAll(keys);

        ///
        /// ↓更新処理ここから
        ///
      

        Matrix4x4 cameraMatrix = MakeAffineMatrix({1.0f, 1.0f, 1.0f}, cameraRotate, cameraTranslate);
        Matrix4x4 viewMatrix = Inverse(cameraMatrix);
        Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
        Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
        Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f); 
        if (isStarted) {
            Vector3 diff = ball.position - spring.anchor;
            float length = Length(diff);
            if (length != 0.0f) {
                Vector3 direction = Normalize(diff);
                Vector3 restPosition = spring.anchor + direction * spring.naturalLength;
                Vector3 displacement = ball.position - restPosition;
                Vector3 restoringForce = -spring.stiffness * displacement;

                // 
                Vector3 dampingForce = -spring.dampingCoefficient * ball.velocity;

                //
                Vector3 force = restoringForce + dampingForce;
                ball.acceleration = force / ball.mass;
            }

            // 加速度と速度を用いてボールの位置と速度を更新する
            ball.velocity += ball.acceleration * deltaTime;
            ball.position += ball.velocity * deltaTime;
        }
        ///
        /// ↑更新処理ここまで
        ///

        ///
        /// ↓描画処理ここから
        ///
       DrawSphere({ball.position, ball.radius}, viewProjectionMatrix, viewportMatrix, ball.color);
        DrawSegment({spring.anchor, ball.position - spring.anchor}, viewProjectionMatrix, viewportMatrix, 0xFFFFFFFF);

        DrawGrid(viewProjectionMatrix, viewportMatrix);
        ImGui::Begin("Window");
     
           if (ImGui::Button("Start")) {
              Reset();
            isStarted = true; 
           
        }
        ImGui::End();
        ///
        /// ↑描画処理ここまで
        ///

        // フレームの終了
        Novice::EndFrame();

        // ESCキーが押されたらループを抜ける
        if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
            break;
        }
    }

    // ライブラリの終了
    Novice::Finalize();
    return 0;
}