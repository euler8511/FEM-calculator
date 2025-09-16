# --- Gmsh 생성 클래스 ---
import sys
import gmsh
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QVBoxLayout, QLabel, QMessageBox

# --- Gmsh 생성 클래스 ---

class GmshGenerator:
    """사용자 입력을 기반으로 Gmsh 메쉬를 생성, 관리 및 내보냅니다."""
    def __init__(self, system_data, force_data, fix_data, youngs_modulus, poisson_ratio):
        self.system_data = system_data
        self.force_data = force_data
        self.fix_data = fix_data
        self.youngs_modulus = youngs_modulus
        self.poisson_ratio = poisson_ratio
        self.output_filename = "generated_mesh.msh"

    def generate_mesh(self):
        """제공된 데이터를 사용하여 메쉬 생성 프로세스를 실행합니다."""
        try:
            gmsh.initialize()
            gmsh.model.add("fem_model")

            # 1. 시스템 데이터로 Box 생성
            x_dim = self.system_data.get('x', 0)
            y_dim = self.system_data.get('y', 0)
            z_dim = self.system_data.get('z', 0)
            mesh_size = self.system_data.get('mesh', 1.0)

            gmsh.model.occ.addBox(0, 0, 0, x_dim, y_dim, z_dim)
            gmsh.model.occ.synchronize()

            # 2. 힘과 고정 위치에 포인트(노드) 생성 및 Physical Group 준비
            tool_points = []
            force_point_tags = []
            fix_point_tags = []

            # 힘 위치 포인트 추가
            for item in self.force_data:
                if not isinstance(item, dict): continue
                px, py, pz = item.get('force_x_pstn'), item.get('force_y_pstn'), item.get('force_z_pstn')
                if all(v is not None for v in [px, py, pz]):
                    pt_tag = gmsh.model.occ.addPoint(px, py, pz)
                    tool_points.append((0, pt_tag))
                    force_point_tags.append(pt_tag)

            # 고정 위치 포인트 추가
            for item in self.fix_data:
                if not isinstance(item, dict): continue
                px, py, pz = item.get('pos_x'), item.get('pos_y'), item.get('pos_z')
                if all(v is not None for v in [px, py, pz]):
                    pt_tag = gmsh.model.occ.addPoint(px, py, pz)
                    tool_points.append((0, pt_tag))
                    fix_point_tags.append(pt_tag)
            
            # 3. Fragment를 사용하여 Box에 포인트 각인
            if tool_points:
                volumes = gmsh.model.getEntities(3)
                gmsh.model.occ.fragment(volumes, tool_points)
                gmsh.model.occ.synchronize()
                print(f"Embedded {len(tool_points)} points into the volume.")

            # 4. Physical Groups 생성
            volumes = gmsh.model.getEntities(3)
            if volumes:
                volume_tags = [v[1] for v in volumes]
                gmsh.model.addPhysicalGroup(3, volume_tags, name="box")
            if force_point_tags:
                gmsh.model.addPhysicalGroup(0, force_point_tags, name="Neumann_BCs")
            if fix_point_tags:
                gmsh.model.addPhysicalGroup(0, fix_point_tags, name="Diri_BCs")

            # 5. 메쉬 옵션 설정 및 생성
            gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
            gmsh.option.setNumber("Mesh.Algorithm3D", 1)
            
            # === [ 중요 수정 사항 ] ===
            # 기존 gmsh.model.mesh.setOrder(2) 대신 전역 옵션을 사용하여 2차 요소를 강제합니다.
            gmsh.option.setNumber("Mesh.ElementOrder", 2)
            
            # (선택) 불완전한 2차 요소가 생성되는 것을 방지하는 옵션
            gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)

            # 3차원 메쉬 생성
            gmsh.model.mesh.generate(3)

            # 6. 파일 저장 및 시각화
            gmsh.write(self.output_filename)
            print(f"Successfully generated '{self.output_filename}'")
            if '-nopopup' not in sys.argv:
                gmsh.option.setNumber("Geometry.Points", 1)
                gmsh.fltk.run()
            
            return True

        except Exception as e:
            error_message = f"An error occurred during Gmsh generation:\n{e}"
            print(error_message)
            QMessageBox.critical(None, "Gmsh Error", error_message)
            return False
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()