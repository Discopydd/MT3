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

    const float moveSpeed = 0.01f;

    Vector3 cameraTranslate = {0.0f, 2.64f, -9.27f};
    Vector3 cameraRotate = {0.26f, 0.0f, 0.0f};
    Vector3 viewTranslate = {0.0f, 0.0f, 0.0f};
    Vector3 cameraScale = {1.0f, 1.0f, 1.0f};

     AABB aabb1 = { {-0.5f, -0.5f, -0.5f}, {0.0f, 0.0f, 0.0f} };
       Sphere sphere ={ {0.0f, 0.0f, 0.0f}, 0.5f };

    bool isDraggingMiddle = false;
    int lastMousePosX = 0, lastMousePosY = 0;
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
        ///カメラ移動
        if (keys[DIK_W]) {
            cameraTranslate.z += moveSpeed;
        }
        if (keys[DIK_S]) {
            cameraTranslate.z -= moveSpeed;
        }
        if (keys[DIK_A]) {
            cameraTranslate.x -= moveSpeed;
        }
        if (keys[DIK_D]) {
            cameraTranslate.x += moveSpeed;
        }
        if (keys[DIK_UP]) {
            cameraTranslate.y += moveSpeed;
        }
        if (keys[DIK_DOWN]) { cameraTranslate.y -= moveSpeed;}
        int mouseX, mouseY;
        Novice::GetMousePosition(&mouseX, &mouseY);
        if (Novice::IsPressMouse(2)) { 
            if (!isDraggingMiddle) {
                isDraggingMiddle = true;
                lastMousePosX = mouseX;
                lastMousePosY = mouseY;
            } else {
                int deltaX = mouseX - lastMousePosX;
                int deltaY = mouseY - lastMousePosY;

                cameraRotate.y += deltaX * 0.001f;
                cameraRotate.x += deltaY * 0.001f;

                lastMousePosX = mouseX;
                lastMousePosY = mouseY;
            }
        } else {
            isDraggingMiddle = false;
        }

        int wheel = Novice::GetWheel();
        Vector3 forward = {
            -sin(cameraRotate.y) * cos(cameraRotate.x),
            sin(cameraRotate.x),
            -cos(cameraRotate.y) * cos(cameraRotate.x)
        };

        cameraTranslate.x -= forward.x * wheel * 0.001f;
        cameraTranslate.y -= forward.y * wheel * 0.001f;
        cameraTranslate.z -= forward.z * wheel * 0.001f;
        
         Matrix4x4 viewMatrix = MakeAffineMatrix({1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f}, viewTranslate);
        Matrix4x4 cameraMatrix = MakeAffineMatrix(cameraScale, cameraRotate, cameraTranslate);
        Matrix4x4 combinedViewMatrix = Multiply(viewMatrix, Inverse(cameraMatrix));
        Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
        Matrix4x4 viewProjectionMatrix = Multiply(combinedViewMatrix, projectionMatrix);
        Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

        bool collision = IsCollisionBaBo(aabb1, sphere);
        uint32_t color1 = collision ? RED : WHITE;
        uint32_t color2 = collision ? RED : WHITE;;
        ///
        /// ↑更新処理ここまで
        ///

        ///
        /// ↓描画処理ここから
        ///

        DrawAABB(aabb1, viewProjectionMatrix, viewportMatrix, color1);
         DrawSphere(sphere, viewProjectionMatrix, viewportMatrix, color2);


        DrawGrid(viewProjectionMatrix, viewportMatrix);
        ImGui::Begin("Window");
        ImGui::DragFloat3("AABB Min", &aabb1.min.x, 0.01f);
        ImGui::DragFloat3("AABB Max", &aabb1.max.x, 0.01f);
        ImGui::DragFloat3("Sphere Center", &sphere.center.x, 0.01f);
        ImGui::DragFloat("Sphere Radius", &sphere.radius, 0.01f);
         aabb1.min.x = min(aabb1.min.x, aabb1.max.x);
        aabb1.max.x = max(aabb1.min.x, aabb1.max.x);
        aabb1.min.y = min(aabb1.min.y, aabb1.max.y);
        aabb1.max.y = max(aabb1.min.y, aabb1.max.y);
        aabb1.min.z = min(aabb1.min.z, aabb1.max.z);
        aabb1.max.z = max(aabb1.min.z, aabb1.max.z);
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