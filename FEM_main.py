# FEM_main.py

import sys
import gmsh
from PyQt5 import uic
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QVBoxLayout, QLabel, QMessageBox
from PyQt5.QtCore import Qt, QStringListModel
from gmsh_creation import GmshGenerator
# 'ReactionSolver'를 Canvas의 파일 이름인 'force_analysis'로 수정했습니다.
from ReactionSolver import ForceAnalysis
import tkinter as tk
from tkinter import filedialog
import meshio
from BeamSolver import BeamAnalysisWindow


# --- UI File Loading ---
# 스크립트와 동일한 디렉토리에 .ui 파일이 있는지 확인합니다.
try:
    Ui_Dialog, QtBaseClass = uic.loadUiType('FEM_calc.ui')
    Ui_Dialog1, QtBaseClass1 = uic.loadUiType('reaction_force.ui')
    Ui_Dialog2, QtBaseClass2 = uic.loadUiType('Beam_analysis.ui')
    Ui_Dialog3, QtBaseClass3 = uic.loadUiType('modal.ui')
    
except FileNotFoundError as e:
    # 필수 UI 파일이 없을 경우, 사용자에게 알리고 프로그램을 종료합니다.
    # 이 부분은 GUI가 시작되기 전이므로 QMessageBox 대신 print를 사용합니다.
    print(f"Error: Could not find a required UI file: {e.filename}")
    print("Please make sure FEM_calc.ui, reaction_force.ui, static.ui, and modal.ui are present.")
    sys.exit(1)

# --- 데이터 입력을 위한 대화상자 클래스 ---

class SystemInfoDialog(QDialog):
    """system_dialogue.ui를 로드하여 시스템 크기 및 메시 입력을 받습니다."""
    def __init__(self, parent=None, initial_data=None):
        super().__init__(parent)
        try:
            uic.loadUi('system_dialogue.ui', self)
        except FileNotFoundError:
            # UI 파일이 없을 경우를 대비한 폴백 UI
            self.setWindowTitle("System Information")
            self.layout = QVBoxLayout()
            self.layout.addWidget(QLabel("Could not load system_dialogue.ui"))
            self.setLayout(self.layout)
            # UI 로딩 실패 시 더 이상 진행하지 않도록 return
            return

        if initial_data:
            self.lineEdit.setText(str(initial_data.get('x', '')))
            self.lineEdit_2.setText(str(initial_data.get('y', '')))
            self.lineEdit_3.setText(str(initial_data.get('z', '')))
            self.lineEdit_4.setText(str(initial_data.get('mesh', '')))


class ForceDialog(QDialog):
    """force.ui를 로드하여 힘 및 위치 입력을 받습니다."""
    def __init__(self, parent=None, initial_data=None):
        super().__init__(parent)
        try:
            uic.loadUi('force.ui', self)
        except FileNotFoundError:
            self.setWindowTitle("Force Information")
            self.layout = QVBoxLayout()
            self.layout.addWidget(QLabel("Could not load force.ui"))
            self.setLayout(self.layout)
            return
            
        if initial_data:
            self.force_x.setText(str(initial_data.get('force_x', '')))
            self.force_y.setText(str(initial_data.get('force_y', '')))
            self.force_z.setText(str(initial_data.get('force_z', '')))
            self.force_x_pstn.setText(str(initial_data.get('force_x_pstn', '')))
            self.force_y_pstn.setText(str(initial_data.get('force_y_pstn', '')))
            self.force_z_pstn.setText(str(initial_data.get('force_z_pstn', '')))


# [NEW] fix.ui를 위한 대화상자 클래스
class FixDialog(QDialog):
    """fix.ui를 로드하여 고정 위치 및 자유도를 지정합니다."""
    def __init__(self, parent=None, initial_data=None):
        super().__init__(parent)
        try:
            uic.loadUi('fix.ui', self)
        except FileNotFoundError:
            self.setWindowTitle("Fix Information")
            self.layout = QVBoxLayout()
            self.layout.addWidget(QLabel("Could not load fix.ui"))
            self.setLayout(self.layout)
            return

        if initial_data:
            # 기존 데이터 편집 시 필드 채우기
            self.fix_x.setText(str(initial_data.get('pos_x', '')))
            self.fix_y.setText(str(initial_data.get('pos_y', '')))
            self.fix_z.setText(str(initial_data.get('pos_z', '')))
            # 값이 0(고정)인지 여부에 따라 체크박스 설정
            self.checkBox_x.setChecked(initial_data.get('fix_x') == 0)
            self.checkBox_y.setChecked(initial_data.get('fix_y') == 0)
            self.checkBox_z.setChecked(initial_data.get('fix_z') == 0)


# --- 메인 애플리케이션 창 ---

class ReactionForceCalculatorWindow(QDialog, Ui_Dialog1):
    """'Reaction Force Calculator' 메인 창."""
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        self.list_model = QStringListModel()
        self.listView.setModel(self.list_model)

        # --- [MODIFIED] 이미지에 명시된 데이터로 사전 초기화 ---
        self.system_data = {
            'x': 0.8, 'y': 0.2, 'z': 0.8, 'mesh': 0.05
        }
        self.force_data_list = [
            {'force_x': 0.0, 'force_y': 3000.0, 'force_z': 0.0,
             'force_x_pstn': 0.4, 'force_y_pstn': 0.2, 'force_z_pstn': 0.4}
        ]
        self.fix_data_list = [
            {'pos_x': 0.0, 'pos_y': 0.0, 'pos_z': 0.0, 'fix_x': 0, 'fix_y': 0, 'fix_z': 0},
            {'pos_x': 0.0, 'pos_y': 0.0, 'pos_z': 0.8, 'fix_x': 0, 'fix_y': 0, 'fix_z': 0},
            {'pos_x': 0.8, 'pos_y': 0.0, 'pos_z': 0.0, 'fix_x': 0, 'fix_y': 0, 'fix_z': 0},
            {'pos_x': 0.8, 'pos_y': 0.0, 'pos_z': 0.8, 'fix_x': 0, 'fix_y': 0, 'fix_z': 0},
        ]
        self.youngs_modul = 2e11
        self.poisson_ratio = 0.3
        self.mesh_file = "generated_mesh.msh"
        self.analysis_instance = None

        # UI 필드에 사전 입력된 값 표시
        self.young_input.setText(str(self.youngs_modul))
        self.poisson_input.setText(str(self.poisson_ratio))
        
        # 리스트 뷰 새로고침
        self._refresh_list_view()
        # ---------------------------------------------------------

        # 시그널과 슬롯 연결
        self.system_info_button.clicked.connect(self.show_system_info_dialog)
        self.force_add_button.clicked.connect(self.show_force_add_dialog)
        self.fix_add_button.clicked.connect(self.show_fix_add_dialog)
        self.edit_button.clicked.connect(self.show_edit_dialog)
        self.remove_button.clicked.connect(self.remove_selected_item)
        self.mesh_update_button.clicked.connect(self.show_mesh_update_dialog)
        self.run_button.clicked.connect(self.run_analysis)
        self.plot_button.clicked.connect(self.plot_results)
        
        
    def _refresh_list_view(self):
        """현재 모든 데이터를 기반으로 QListView를 새로고침하는 중앙 함수."""
        info_items = []
        # 1. 시스템 정보가 있으면 추가
        if self.system_data:
            sys_size_str = f"System Size: ({self.system_data['x']}, {self.system_data['y']}, {self.system_data['z']}), Mesh Size: {self.system_data['mesh']}"
            info_items.append(sys_size_str)

        # 2. 모든 힘 정보 추가
        for i, force_data in enumerate(self.force_data_list):
            force_str = (f"Force {i+1}: ({force_data['force_x']}, {force_data['force_y']}, {force_data['force_z']}) "
                         f"@ Pos: ({force_data['force_x_pstn']}, {force_data['force_y_pstn']}, {force_data['force_z_pstn']})")
            info_items.append(force_str)
            
        # 3. [NEW] 모든 고정 정보 추가
        for i, fix_data in enumerate(self.fix_data_list):
            fixed_axes = []
            if fix_data['fix_x'] == 0: fixed_axes.append('X')
            if fix_data['fix_y'] == 0: fixed_axes.append('Y')
            if fix_data['fix_z'] == 0: fixed_axes.append('Z')
            axes_str = ', '.join(fixed_axes) if fixed_axes else "None"

            fix_str = (f"Fix {i+1}: Pos ({fix_data['pos_x']}, {fix_data['pos_y']}, {fix_data['pos_z']}) "
                       f"- Fixed Dof: [{axes_str}]")
            info_items.append(fix_str)

        # 4. 모델을 업데이트하여 UI 새로고침
        self.list_model.setStringList(info_items)

    def show_system_info_dialog(self):
        """시스템 정보 추가 또는 편집을 위한 대화상자 표시 (업데이트)."""
        dialog = SystemInfoDialog(self, initial_data=self.system_data)
        if dialog.exec_() == QDialog.Accepted:
            try:
                # [ADDED] 사용자가 숫자가 아닌 값을 입력할 경우를 대비한 예외 처리
                self.system_data = {
                    'x': float(dialog.lineEdit.text()),
                    'y': float(dialog.lineEdit_2.text()),
                    'z': float(dialog.lineEdit_3.text()),
                    'mesh': float(dialog.lineEdit_4.text())
                }
                print("--- System Information Updated ---")
                self._refresh_list_view()
            except ValueError:
                QMessageBox.warning(self, "Input Error", "잘못된 입력입니다. 모든 필드에 유효한 숫자를 입력하세요.")
                print("Error: Invalid numeric input for system information.")
        else:
            print("System information input canceled.")

    def show_force_add_dialog(self):
        """새로운 힘 정보를 추가하기 위한 대화상자 표시 (추가)."""
        dialog = ForceDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            try:
                # [ADDED] 숫자가 아닌 값 입력에 대한 예외 처리
                new_force = {
                    'force_x': float(dialog.force_x.text()),
                    'force_y': float(dialog.force_y.text()),
                    'force_z': float(dialog.force_z.text()),
                    'force_x_pstn': float(dialog.force_x_pstn.text()),
                    'force_y_pstn': float(dialog.force_y_pstn.text()),
                    'force_z_pstn': float(dialog.force_z_pstn.text())
                }
                self.force_data_list.append(new_force)
                print("--- Force Information Added ---")
                self._refresh_list_view()
            except ValueError:
                QMessageBox.warning(self, "Input Error", "잘못된 입력입니다. 모든 힘 관련 필드에 유효한 숫자를 입력하세요.")
                print("Error: Invalid numeric input for force information.")
        else:
            print("Force input canceled.")
    
    # [NEW] 'Fix Add' 버튼을 위한 구현된 메소드
    def show_fix_add_dialog(self):
        """새로운 고정 정보를 추가하기 위한 대화상자 표시 (추가)."""
        dialog = FixDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            try:
                # [ADDED] 숫자가 아닌 값 입력에 대한 예외 처리
                new_fix = {
                    'pos_x': float(dialog.fix_x.text()),
                    'pos_y': float(dialog.fix_y.text()),
                    'pos_z': float(dialog.fix_z.text()),
                    # 체크박스가 선택되면 값은 0(고정), 아니면 None
                    'fix_x': 0 if dialog.checkBox_x.isChecked() else None,
                    'fix_y': 0 if dialog.checkBox_y.isChecked() else None,
                    'fix_z': 0 if dialog.checkBox_z.isChecked() else None,
                }
                self.fix_data_list.append(new_fix)
                print("--- Fix Information Added ---")
                self._refresh_list_view()
            except ValueError:
                QMessageBox.warning(self, "Input Error", "잘못된 입력입니다. 모든 위치 필드에 유효한 숫자를 입력하세요.")
                print("Error: Invalid numeric input for fix information.")
        else:
            print("Fix input canceled.")

    def show_edit_dialog(self):
        """리스트 뷰에서 선택된 항목을 편집합니다."""
        selected_indexes = self.listView.selectedIndexes()
        if not selected_indexes:
            QMessageBox.information(self, "Selection Required", "편집할 항목을 선택하세요.")
            return

        row = selected_indexes[0].row()
        selected_text = selected_indexes[0].data()

        if "System Size" in selected_text:
            self.show_system_info_dialog()

        elif "Force" in selected_text:
            force_list_index = row - (1 if self.system_data else 0)
            if 0 <= force_list_index < len(self.force_data_list):
                dialog = ForceDialog(self, initial_data=self.force_data_list[force_list_index])
                if dialog.exec_() == QDialog.Accepted:
                    try:
                        # [ADDED] 편집 중 숫자가 아닌 값 입력에 대한 예외 처리
                        updated_force = {
                            'force_x': float(dialog.force_x.text()),'force_y': float(dialog.force_y.text()),'force_z': float(dialog.force_z.text()),
                            'force_x_pstn': float(dialog.force_x_pstn.text()),'force_y_pstn': float(dialog.force_y_pstn.text()),'force_z_pstn': float(dialog.force_z_pstn.text())
                        }
                        self.force_data_list[force_list_index] = updated_force
                        print(f"--- Force {force_list_index + 1} Updated ---")
                        self._refresh_list_view()
                    except ValueError:
                        QMessageBox.warning(self, "Input Error", "잘못된 입력입니다. 모든 힘 관련 필드에 유효한 숫자를 입력하세요.")
                        print("Error: Invalid numeric input while editing force.")
                else:
                    print("Force edit canceled.")
        
        elif "Fix" in selected_text: # [NEW] Fix 항목에 대한 편집 로직
            start_of_fix_items = (1 if self.system_data else 0) + len(self.force_data_list)
            fix_list_index = row - start_of_fix_items
            
            if 0 <= fix_list_index < len(self.fix_data_list):
                dialog = FixDialog(self, initial_data=self.fix_data_list[fix_list_index])
                if dialog.exec_() == QDialog.Accepted:
                    try:
                        # [ADDED] 편집 중 숫자가 아닌 값 입력에 대한 예외 처리
                        updated_fix = {
                            'pos_x': float(dialog.fix_x.text()), 'pos_y': float(dialog.fix_y.text()), 'pos_z': float(dialog.fix_z.text()),
                            'fix_x': 0 if dialog.checkBox_x.isChecked() else None,
                            'fix_y': 0 if dialog.checkBox_y.isChecked() else None,
                            'fix_z': 0 if dialog.checkBox_z.isChecked() else None,
                        }
                        self.fix_data_list[fix_list_index] = updated_fix
                        print(f"--- Fix {fix_list_index + 1} Updated ---")
                        self._refresh_list_view()
                    except ValueError:
                        QMessageBox.warning(self, "Input Error", "잘못된 입력입니다. 모든 위치 필드에 유효한 숫자를 입력하세요.")
                        print("Error: Invalid numeric input while editing fix.")
                else:
                    print("Fix edit canceled.")


    def remove_selected_item(self):
        """리스트 뷰와 데이터 저장소에서 선택된 항목을 제거합니다."""
        selected_indexes = self.listView.selectedIndexes()
        if not selected_indexes:
            QMessageBox.information(self, "Selection Required", "삭제할 항목을 선택하세요.")
            return
            
        row = selected_indexes[0].row()
        selected_text = selected_indexes[0].data()

        if "System Size" in selected_text:
            self.system_data.clear()
            print("--- System Information Removed ---")
        elif "Force" in selected_text:
            force_list_index = row - (1 if self.system_data else 0)
            if 0 <= force_list_index < len(self.force_data_list):
                self.force_data_list.pop(force_list_index)
                print(f"--- Force {force_list_index + 1} Removed ---")
        elif "Fix" in selected_text: # [NEW] Fix 항목에 대한 제거 로직
            start_of_fix_items = (1 if self.system_data else 0) + len(self.force_data_list)
            fix_list_index = row - start_of_fix_items
            if 0 <= fix_list_index < len(self.fix_data_list):
                self.fix_data_list.pop(fix_list_index)
                print(f"--- Fix {fix_list_index + 1} Removed ---")

        self._refresh_list_view()

    def show_mesh_update_dialog(self):
        try:
            # [ADDED] 숫자가 아닌 값 입력에 대한 예외 처리
            self.youngs_modul = float(self.young_input.text())
            self.poisson_ratio = float(self.poisson_input.text())
            print("--- Material Properties Updated ---")
            print(f"Young's Modulus: {self.youngs_modul}, Poisson's Ratio: {self.poisson_ratio}")
            print("--- Current Data ---")
            print("System Data:", self.system_data)
            print("Force Data:", self.force_data_list)
            print("Fix Data:", self.fix_data_list)
            print("=====================================")
            self.gmsh_ = GmshGenerator(self.system_data, self.force_data_list, self.fix_data_list, 
                                       self.youngs_modul, self.poisson_ratio)
            self.gmsh_.generate_mesh()
            
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Young's Modulus 또는 Poisson's Ratio에 유효한 숫자를 입력하세요.")
            print("Error: Invalid input for material properties.")


    def run_analysis(self):
        # [MODIFIED] 다중 힘을 처리하도록 수정
        print("--- Starting FEM Analysis ---")
        try:
            # 1. 필수 데이터 확인
            if not all([self.youngs_modul, self.poisson_ratio, self.force_data_list]):
                QMessageBox.warning(self, "Data Missing", "해석을 실행하려면 재료 물성치와 하나 이상의 힘 조건이 필요합니다.")
                return

            # 2. ForceAnalysis 인스턴스 생성 및 실행
            # force_data_list 전체를 전달합니다.
            self.analysis_instance = ForceAnalysis(
                msh_file=self.mesh_file,
                force_data=self.force_data_list,  # [CHANGED]
                fix_data=self.fix_data_list,
                E=self.youngs_modul,
                v=self.poisson_ratio
            )
            self.analysis_instance.run_simulation()

            QMessageBox.information(self, "Success", "해석이 완료되었습니다. Plot 버튼으로 결과를 확인하세요.")
            print("--- FEM Analysis Finished Successfully ---")

        except Exception as e:
            error_message = f"해석 중 오류가 발생했습니다:\n{e}"
            print(error_message)
            QMessageBox.critical(self, "Analysis Error", error_message)
            self.analysis_instance = None


    def plot_results(self):
        # [NEW] This function plots the results from the completed analysis
        if self.analysis_instance is None:
            QMessageBox.warning(self, "No Results", "먼저 해석을 실행해야 합니다 (RUN 버튼).")
            return
        
        try:
            self.analysis_instance.plot()
        except Exception as e:
            error_message = f"결과 시각화 중 오류가 발생했습니다:\n{e}"
            print(error_message)
            QMessageBox.critical(self, "Plotting Error", error_message)


class CustomDialog(QDialog):
    """간단한 메시지를 표시하기 위한 일반적인 대화상자."""
    def __init__(self, title, message, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setGeometry(400, 400, 300, 100)
        layout = QVBoxLayout()
        label = QLabel(message)
        layout.addWidget(label)
        self.setLayout(layout)
        
                
        
class ShaftAnalysisWindow(QDialog, Ui_Dialog3):
    """'Modal Analysis' 창 (플레이스홀더)."""
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        
        
class Static2DAnalysisWindow(QDialog, Ui_Dialog3):
    """'Modal Analysis' 창 (플레이스홀더)."""
    def __init__(self):
        super().__init__()
        self.setupUi(self)


class PipeThermalAnalysisWindow(QDialog, Ui_Dialog3):
    """'Modal Analysis' 창 (플레이스홀더)."""
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        
        
        
# --- 메인 선택 대화상자 ---

class SelectionDialog(QDialog, Ui_Dialog):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.new_window = None
        self.buttonBox.accepted.connect(self.on_ok_button_clicked)

    def on_ok_button_clicked(self):
        """콤보박스 선택에 따라 적절한 메인 창을 생성합니다."""
        selected_option = self.comboBox.currentText()
        if selected_option == "Reaction Force Calculator with GMSH":
            self.new_window = ReactionForceCalculatorWindow()
        elif selected_option == "Beam Analysis(Timoshenko beam) with GMSH":
            self.new_window = BeamAnalysisWindow()
        elif selected_option == "Shaft modal with GMSH":
            self.new_window = ShaftAnalysisWindow()
        elif selected_option == "2D Static Analysis with GMSH":
            self.new_window = Static2DAnalysisWindow()
        elif selected_option == "Pipe Thermal Stress Analysis with GMSH":
            self.new_window = PipeThermalAnalysisWindow()

# --- 프로그램 실행 ---

if __name__ == '__main__':
    app = QApplication(sys.argv)
    
    selection_dialog = SelectionDialog()
    
    # 초기 선택 대화상자 표시
    if selection_dialog.exec_() == QDialog.Accepted and selection_dialog.new_window:
        # 사용자가 OK를 클릭하면 선택된 메인 창을 표시
        selection_dialog.new_window.show()
        sys.exit(app.exec_())
    
    # 사용자가 취소하면 종료
    sys.exit(0)