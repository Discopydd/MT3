#include <Novice.h>
#include<imgui.h>
#include"MyMath.h"

static const int kWindowWidth = 1280; 
static const int kWindowHeight = 720;

const char kWindowTitle[] = "LE2C_29_リ_ヨン";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	Vector3 cameraTranslate{0.0f, 1.9f, -6.49f}; 
Vector3 cameraRotate{0.26f, 0.0f, 0.0f};
Segment segment = {{-2.0f, -1.0f, 0.0f}, {3.0f, 2.0f, 2.0f}};
Vector3 point = {-1.5f, 0.6f, 0.6f};
 Sphere sphere = {{0.0f, 0.0f, 0.0f}, 1.0f};

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
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

        Vector3 closestPoint = ClosestPoint(point, segment);

        Sphere pointSphere = {point, 0.01f};
        Sphere closestPointSphere = {closestPoint, 0.01f};

		  Vector3 start = Transform(Transform(segment.origir, viewProjectionMatrix), viewportMatrix);
    Vector3 end = Transform(Transform(Add(segment.origir,segment.diff),viewProjectionMatrix), viewportMatrix);

		Vector3 project = Project(Subtract(point, segment.origir), segment.diff);

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///
		  DrawGrid(viewProjectionMatrix, viewportMatrix);
		    DrawSphere(pointSphere, viewProjectionMatrix, viewportMatrix, 0xFF0000FF);  // Red
        DrawSphere(closestPointSphere, viewProjectionMatrix, viewportMatrix, 0x000000FF);  // Black
		    Novice::DrawLine((int)start.x, (int)start.y,
                     (int)end.x, (int)end.y, 0xFFFFFFFF);
		ImGui::Begin("Window");
        ImGui::InputFloat3("Point", &point.x, "%.3f", ImGuiInputTextFlags_ReadOnly);
        ImGui::InputFloat3("Segment origin", &segment.origir.x, "%.3f", ImGuiInputTextFlags_ReadOnly);
        ImGui::InputFloat3("Segment diff", &segment.diff.x, "%.3f", ImGuiInputTextFlags_ReadOnly);
		ImGui::InputFloat3("Project", &project.x, "%.3f", ImGuiInputTextFlags_ReadOnly);
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
