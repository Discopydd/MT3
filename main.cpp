#include <Novice.h>
#include<imgui.h>
#include"MyMath.h"
#include "Vector3.h"
#include "Matrix4x4.h"


static const int kWindowWidth = 1280; 
static const int kWindowHeight = 720;

const char kWindowTitle[] = "LE2C_29_リ_ヨン";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

    // ライブラリの初期化
    Novice::Initialize(kWindowTitle, 1280, 720);

 

    Vector3 cameraTranslate = {0.0f, 2.64f, -9.27f};
    Vector3 cameraRotate = {0.26f, 0.0f, 0.0f};
    Vector3 viewTranslate = {0.0f, 0.0f, 0.0f};
    Vector3 cameraScale = {1.0f, 1.0f, 1.0f};

  
    Vector3 translates[3] = {
        {0.2f, 1.0f, 0.0f}, // 肩
        {0.4f, 0.0f, 0.0f}, // 肘
        {0.3f, 0.0f, 0.0f}  // 手
    };

    Vector3 rotates[3] = {
        {0.0f, 0.0f, -6.8f}, // 肩
        {0.0f, 0.0f, -1.4f}, // 肘
        {0.0f, 0.0f, 0.0f}   // 手
    };

    Vector3 scales[3] = {
        {1.0f, 1.0f, 1.0f}, // 肩
        {1.0f, 1.0f, 1.0f}, // 肘
        {1.0f, 1.0f, 1.0f}  // 手
    };

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
        
         Matrix4x4 shoulderLocalMatrix = MakeAffineMatrix(scales[0], rotates[0], translates[0]);
        Matrix4x4 shoulderWorldMatrix = shoulderLocalMatrix;

        Matrix4x4 elbowLocalMatrix = MakeAffineMatrix(scales[1], rotates[1], translates[1]);
        Matrix4x4 elbowWorldMatrix = Multiply(elbowLocalMatrix, shoulderWorldMatrix);

        Matrix4x4 handLocalMatrix = MakeAffineMatrix(scales[2], rotates[2], translates[2]);
        Matrix4x4 handWorldMatrix = Multiply(handLocalMatrix, elbowWorldMatrix);


        Sphere shoulder = { Transform({0.0f, 0.0f, 0.0f}, shoulderWorldMatrix), 0.05f };
        Sphere elbow = { Transform({0.0f, 0.0f, 0.0f}, elbowWorldMatrix), 0.05f };
        Sphere hand = { Transform({0.0f, 0.0f, 0.0f}, handWorldMatrix), 0.05f };

       
        ///
        /// ↑更新処理ここまで
        ///

        ///
        /// ↓描画処理ここから
        ///
       DrawSphere(shoulder, viewProjectionMatrix, viewportMatrix, RED);
        DrawSphere(elbow, viewProjectionMatrix, viewportMatrix, GREEN);
        DrawSphere(hand, viewProjectionMatrix, viewportMatrix, BLUE);
        // 肩-肘、肘-手の線を描画
        Segment shoulderToElbow = {shoulder.center, Subtract(elbow.center, shoulder.center)};
        Segment elbowToHand = {elbow.center, Subtract(hand.center, elbow.center)};

        DrawSegment(shoulderToElbow, viewProjectionMatrix, viewportMatrix, WHITE);
        DrawSegment(elbowToHand, viewProjectionMatrix, viewportMatrix, WHITE);
        DrawGrid(viewProjectionMatrix, viewportMatrix);
        ImGui::Begin("Window");
      ImGui::DragFloat3("translates[0]", &translates[0].x, 0.01f);
         ImGui::DragFloat3("rotates[0]", &rotates[0].x, 0.01f);
           ImGui::DragFloat3("scales[0]", &scales[0].x, 0.01f);
        ImGui::DragFloat3("translates[1]", &translates[1].x, 0.01f);
         ImGui::DragFloat3("rotates[1]", &rotates[1].x, 0.01f);
           ImGui::DragFloat3("scales[1]", &scales[1].x, 0.01f);
        ImGui::DragFloat3("translates[2]", &translates[2].x, 0.01f);
        ImGui::DragFloat3("rotates[2]", &rotates[2].x, 0.01f);
        ImGui::DragFloat3("scales[2]", &scales[2].x, 0.01f);
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