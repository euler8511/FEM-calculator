import sys
import numpy as np
import matplotlib.pyplot as plt
import meshio
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
from tkinter import filedialog

from PyQt5 import uic
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QVBoxLayout, QLabel, QMessageBox
from PyQt5.QtCore import QStringListModel
from PyQt5.QtGui import QPixmap

# Import sectionproperties library
import sectionproperties.pre.library as lib
from sectionproperties.analysis import Section
import section_rc

# --- UI File Loading ---
try:
    Ui_Dialog, QtBaseClass = uic.loadUiType('Beam_analysis.ui')
except FileNotFoundError as e:
    print(f"Error: Could not find a required UI file: {e.filename}")
    print("Please make sure all .ui files are present in the same directory.")
    sys.exit(1)

# --- MODIFIED FUNCTION for calculating section properties using sectionproperties ---
def calculate_section_properties(section_type, params):
    """
    Calculates cross-sectional properties using the sectionproperties library.
    params: A dictionary of parameters for the section.
    Returns: A, I_y, I_z, J, kappa_y, kappa_z, c_y_max, c_z_max
    """
    geometry = None
    try:
        # Create geometry based on section type
        if section_type == "hollow box section":
            geometry = lib.rectangular_hollow_section(
                d=params['d'], b=params['b'], t=params['t'], r_out=params.get('r_out', 0), n_r=8
            )
        elif section_type == "I section":
            geometry = lib.i_section(
                d=params['d'], b=params['b'], t_f=params['t_f'], t_w=params['t_w'], r=params.get('r', 0), n_r=8
            )
        elif section_type == "C section":
            geometry = lib.channel_section(
                d=params['d'], b=params['b'], t_f=params['t_f'], t_w=params['t_w'], r=params.get('r', 0), n_r=8
            )
        elif section_type == "L section":
            geometry = lib.angle_section(
                d=params['d'], b=params['b'], t=params['t'], r_r=params.get('r_r', 0), r_t=params.get('r_t', 0), n_r=8
            )
        elif section_type == "rectangular section":
            geometry = lib.rectangular_section(d=params['d'], b=params['b'])
        elif section_type == "circular section":
            geometry = lib.circular_section(d=params['d'], n=64)
        elif section_type == "hollow circular section":
            geometry = lib.circular_hollow_section(d=params['d'], t=params['t'], n=64)
        else:
            print(f"Warning: Unknown section type '{section_type}'.")
            return (0,) * 8

        # Create a mesh and Section object.
        t_vals = [v for k, v in params.items() if 't' in k and v > 0]
        mesh_size_ref = min(t_vals) if t_vals else params.get('d', 1)
        mesh_size = mesh_size_ref / 4.0 if mesh_size_ref > 0 else 0.1
        
        geometry.create_mesh(mesh_sizes=[mesh_size])
        section = Section(geometry, time_info=False)

        # Perform analysis
        section.calculate_geometric_properties()
        section.calculate_warping_properties()

        # Extract properties
        A = section.get_area()
        I_x, I_y, I_xy = section.get_ic()
        J = section.get_j()

        shear_area_x, shear_area_y = section.get_as()
        kappa_y = shear_area_x / A
        kappa_z = shear_area_y / A

        # Compute max distances from centroid
        cx, cy = section.get_c()

        vertices = section.mesh['vertices']
        c_y_max = max(abs(vertices[:, 0] - cx))
        c_z_max = max(abs(vertices[:, 1] - cy))

        
        return A, I_x, I_y, J, kappa_y, kappa_z, c_y_max, c_z_max

    except Exception as e:
        print(f"Error in sectionproperties for {section_type} with params {params}: {e}")
        return (0,) * 8

# --- NEW Section Dialog Classes ---

class BaseSectionDialog(QDialog):
    """Base class for section dialogs to handle common functionality."""
    def __init__(self, ui_file, parent=None, initial_data=None):
        super().__init__(parent)
        try:
            uic.loadUi(ui_file, self)
            self.setWindowTitle(f"Section Information")
            if initial_data:
                self.set_data(initial_data)
        except FileNotFoundError:
            self.setWindowTitle("UI File Error")
            self.setLayout(QVBoxLayout())
            self.layout().addWidget(QLabel(f"Could not load {ui_file}"))

    def get_values(self):
        """Placeholder for getting values. To be implemented by subclasses."""
        raise NotImplementedError

    def set_data(self, data):
        """Placeholder for setting data. To be implemented by subclasses."""
        pass

    def _get_float(self, widget_name, default=0.0):
        """Helper to safely get a float value from a widget."""
        if hasattr(self, widget_name):
            try:
                print("dddddddddddd")
                print(getattr(self, widget_name).text())
                return float(getattr(self, widget_name).text().strip())
            except (ValueError, AttributeError):
                return default
        return default
    
    def _set_text(self, widget_name, value):
        """Helper to safely set text in a widget."""
        if hasattr(self, widget_name) and value is not None:
            getattr(self, widget_name).setText(str(value))

class ISectionDialog(BaseSectionDialog):
    def __init__(self, parent=None, initial_data=None):
        super().__init__('I_section.ui', parent, initial_data)

    def get_values(self):
        return {
            'd': self._get_float('d_input'),
            'b': self._get_float('b_input'),
            't_w': self._get_float('tw_input'),
            't_f': self._get_float('tf_input'),
            'r': self._get_float('r_input'), # Assuming r is not in this UI
        }

    def set_data(self, data):
        self._set_text('d_input', data.get('d'))
        self._set_text('b_input', data.get('b'))
        self._set_text('tw_input', data.get('t_w'))
        self._set_text('tf_input', data.get('t_f'))
        self._set_text('r_input', data.get('r'))

class CSectionDialog(BaseSectionDialog):
    def __init__(self, parent=None, initial_data=None):
        super().__init__('C_section.ui', parent, initial_data)

    def get_values(self):
        return {
            'd': self._get_float('d_input'),
            'b': self._get_float('b_input'),
            't_f': self._get_float('tf_input'),
            't_w': self._get_float('tw_input'),
            'r': self._get_float('r_input'),
        }
    
    def set_data(self, data):
        self._set_text('d_input', data.get('d'))
        self._set_text('b_input', data.get('b'))
        self._set_text('tf_input', data.get('t_f'))
        self._set_text('tw_input', data.get('t_w'))
        self._set_text('r_input', data.get('r'))

class LSectionDialog(BaseSectionDialog):
    def __init__(self, parent=None, initial_data=None):
        super().__init__('L_section.ui', parent, initial_data)
        
    def get_values(self):
        return {
            'd': self._get_float('d_input'),
            'b': self._get_float('b_input'),
            't': self._get_float('t_input'),
            'r_r': self._get_float('r_r_input'),
            'r_t': self._get_float('r_t_input'),
        }

    def set_data(self, data):
        self._set_text('d_input', data.get('d'))
        self._set_text('b_input', data.get('b'))
        self._set_text('t_input', data.get('t'))
        self._set_text('r_r_input', data.get('r_r'))
        self._set_text('r_t_input', data.get('r_t'))

class HollowBoxSectionDialog(BaseSectionDialog):
    def __init__(self, parent=None, initial_data=None):
        super().__init__('hollow_box_section.ui', parent, initial_data)

    def get_values(self):
        return {
            'd': self._get_float('d_input'),
            'b': self._get_float('b_input'),
            't': self._get_float('t_input'),
            'r_out': self._get_float('r_input'), # Assuming r_input is r_out
        }

    def set_data(self, data):
        self._set_text('d_input', data.get('d'))
        self._set_text('b_input', data.get('b'))
        self._set_text('t_input', data.get('t'))
        self._set_text('r_input', data.get('r_out'))

class RectangularSectionDialog(BaseSectionDialog):
    def __init__(self, parent=None, initial_data=None):
        super().__init__('Rectangular_section.ui', parent, initial_data)
    
    def get_values(self):
        return {
            'd': self._get_float('d_input'),
            'b': self._get_float('b_input'),
        }

    def set_data(self, data):
        self._set_text('d_input', data.get('d'))
        self._set_text('b_input', data.get('b'))

class CircularSectionDialog(BaseSectionDialog):
    def __init__(self, parent=None, initial_data=None):
        super().__init__('circular_section.ui', parent, initial_data)

    def get_values(self):
        return { 'd': self._get_float('d_input') }

    def set_data(self, data):
        self._set_text('d_input', data.get('d'))

class HollowCircularSectionDialog(BaseSectionDialog):
    def __init__(self, parent=None, initial_data=None):
        super().__init__('hollow_circular_section.ui', parent, initial_data)
    
    def get_values(self):
        return {
            'd': self._get_float('d_input'),
            't': self._get_float('t_input'),
        }

    def set_data(self, data):
        self._set_text('d_input', data.get('d'))
        self._set_text('t_input', data.get('t'))

class BoundaryConditionDialog(QDialog):
    """Loads beam_fix.ui and manages boundary condition data."""
    def __init__(self, parent=None, initial_data=None):
        super().__init__(parent)
        try:
            uic.loadUi('beam_fix.ui', self)
            self.setWindowTitle("Set Boundary Condition")

            if hasattr(self, 'fix_radioButton'):
                self.fix_radioButton.toggled.connect(self._update_widget_states)
            if hasattr(self, 'applyForce_radioButton'):
                self.applyForce_radioButton.toggled.connect(self._update_widget_states)

            if initial_data:
                self.set_data(initial_data)
            else:
                self._update_widget_states()

        except FileNotFoundError:
            self.setWindowTitle("Error")
            self.setLayout(QVBoxLayout())
            self.layout().addWidget(QLabel("Could not load beam_fix.ui"))

    def _update_widget_states(self):
        is_fix = self.fix_radioButton.isChecked()
        is_force = self.applyForce_radioButton.isChecked()

        for name in ['x_checkBox', 'y_checkBox', 'z_checkBox']:
            if hasattr(self, name): getattr(self, name).setEnabled(is_fix)
        for name in ['force_x_input', 'force_y_input', 'force_z_input']:
            if hasattr(self, name): getattr(self, name).setEnabled(is_force)

    def set_data(self, data):
        if data.get('type') == 'Fix':
            self.fix_radioButton.setChecked(True)
            self.x_checkBox.setChecked(data.get('fix_x', False))
            self.y_checkBox.setChecked(data.get('fix_y', False))
            self.z_checkBox.setChecked(data.get('fix_z', False))
        elif data.get('type') == 'Force':
            self.applyForce_radioButton.setChecked(True)
            self.force_x_input.setText(str(data.get('force_x', 0.0)))
            self.force_y_input.setText(str(data.get('force_y', 0.0)))
            self.force_z_input.setText(str(data.get('force_z', 0.0)))
        self._update_widget_states()

    def get_data(self):
        try:
            if self.fix_radioButton.isChecked():
                return {
                    'type': 'Fix',
                    'fix_x': self.x_checkBox.isChecked(),
                    'fix_y': self.y_checkBox.isChecked(),
                    'fix_z': self.z_checkBox.isChecked()
                }
            elif self.applyForce_radioButton.isChecked():
                return {
                    'type': 'Force',
                    'force_x': float(self.force_x_input.text()),
                    'force_y': float(self.force_y_input.text()),
                    'force_z': float(self.force_z_input.text())
                }
            return None
        except ValueError:
            QMessageBox.critical(self, "Input Error", "Please enter valid numbers for forces.")
            return None


class BeamAnalysisWindow(QDialog, Ui_Dialog):
    """Main window for the 3D Timoshenko Beam Analysis application."""
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.setWindowTitle("3D Timoshenko Beam Analysis")

        self.mesh = None
        self.points = None
        self.conn = None
        self.u = None
        self.smoothed_stresses = None
        self.section_data = []
        self.bc_data = []

        self.list_model = QStringListModel()
        self.listView.setModel(self.list_model)

        self.section_type_combo.addItems([
            "I section", "C section", "L section", "hollow box section",
            "rectangular section", "circular section", "hollow circular section"
        ])

        self.mesh_sel_btn.clicked.connect(self.mesh_load)
        self.beam_button.clicked.connect(self.assign_beam_section)
        self.bc_button.clicked.connect(self.assign_bc)
        self.edit_button.clicked.connect(self.edit_item)
        self.remove_button.clicked.connect(self.remove_item)
        self.mesh_update_button.clicked.connect(self.mesh_update)
        self.run_button.clicked.connect(self.run_simulation)
        self.plot_button.clicked.connect(self.plot_results)

    def mesh_load(self):
        root = tk.Tk()
        root.withdraw()
        mesh_path = filedialog.askopenfilename(title="Select Gmsh .msh file", filetypes=[("Gmsh mesh", "*.msh")])
        if not mesh_path: return

        try:
            self.mesh = meshio.read(mesh_path)
            self.root_indicate.setText(mesh_path)
            self.points = self.mesh.points
            self.conn = self.mesh.cells_dict.get('line')
            if self.conn is None: raise ValueError("No 'line' elements in .msh file.")

            physical_groups = list(self.mesh.field_data.keys()) if self.mesh.field_data else []
            self.physical_group_combo.clear()
            self.bc_combo.clear()
            self.physical_group_combo.addItems(physical_groups)
            self.bc_combo.addItems(physical_groups)
            
            self.section_data.clear()
            self.bc_data.clear()
            self.update_list_view()
        except Exception as e:
            QMessageBox.critical(self, "Mesh Load Error", f"Failed to read mesh file:\n{e}")
            self.mesh = None

    def assign_beam_section(self):
        """Assigns section properties to a physical group using dynamic dialogs."""
        group = self.physical_group_combo.currentText()
        if not group:
            QMessageBox.warning(self, "Selection Error", "Please select a physical group.")
            return

        if any(d['group'] == group for d in self.section_data):
            reply = QMessageBox.question(self, 'Confirm Overwrite', f"Overwrite existing section for '{group}'?", QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.No:
                return
            self.section_data = [d for d in self.section_data if d['group'] != group]
        
        section_type = self.section_type_combo.currentText()
        
        dialog_map = {
            "I section": ISectionDialog,
            "C section": CSectionDialog,
            "L section": LSectionDialog,
            "hollow box section": HollowBoxSectionDialog,
            "rectangular section": RectangularSectionDialog,
            "circular section": CircularSectionDialog,
            "hollow circular section": HollowCircularSectionDialog,
        }

        dialog_class = dialog_map.get(section_type)
        if not dialog_class:
            QMessageBox.critical(self, "Error", f"No dialog found for section type '{section_type}'.")
            return
            
        dialog = dialog_class(self)
        if dialog.exec_() == QDialog.Accepted:
            params = dialog.get_values()
            if params is not None:
                self.section_data.append({
                    'group': group, 
                    'type': section_type, 
                    'params': params
                })
                self.update_list_view()

    def assign_bc(self):
        """Assigns a boundary condition to a physical group."""
        group = self.bc_combo.currentText()
        if not group:
            QMessageBox.warning(self, "Selection Error", "Please select a physical group for the BC.")
            return
            
        if any(d['group'] == group for d in self.bc_data):
            reply = QMessageBox.question(self, 'Confirm Overwrite', f"Overwrite existing BC for '{group}'?", QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.No: return
            self.bc_data = [d for d in self.bc_data if d['group'] != group]

        dialog = BoundaryConditionDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            data = dialog.get_data()
            if data:
                data['group'] = group
                self.bc_data.append(data)
                self.update_list_view()

    def edit_item(self):
        """Edits the selected item in the list view."""
        selected = self.listView.selectedIndexes()
        if not selected:
            QMessageBox.warning(self, "Selection Error", "Please select an item to edit.")
            return
        
        row = selected[0].row()
        num_sections = len(self.section_data)

        if row < num_sections: # It's a section
            data = self.section_data[row]
            section_type = data['type']
            initial_params = data['params']

            dialog_map = {
                "I section": ISectionDialog,
                "C section": CSectionDialog,
                "L section": LSectionDialog,
                "hollow box section": HollowBoxSectionDialog,
                "rectangular section": RectangularSectionDialog,
                "circular section": CircularSectionDialog,
                "hollow circular section": HollowCircularSectionDialog,
            }
            dialog_class = dialog_map.get(section_type)
            if not dialog_class: return

            dialog = dialog_class(self, initial_data=initial_params)
            if dialog.exec_() == QDialog.Accepted:
                new_params = dialog.get_values()
                if new_params:
                    data['params'] = new_params
        
        else: # It's a boundary condition
            bc_index = row - num_sections
            data = self.bc_data[bc_index]
            dialog = BoundaryConditionDialog(self, initial_data=data)
            if dialog.exec_() == QDialog.Accepted:
                new_data = dialog.get_data()
                if new_data:
                    new_data['group'] = data['group']
                    self.bc_data[bc_index] = new_data
        
        self.update_list_view()

    def remove_item(self):
        """Removes the selected item from the list view."""
        selected = self.listView.selectedIndexes()
        if not selected:
            QMessageBox.warning(self, "Selection Error", "Please select an item to remove.")
            return

        row = selected[0].row()
        num_sections = len(self.section_data)
        
        if row < num_sections:
            item_text = self.section_data[row]['group']
        else:
            item_text = self.bc_data[row - num_sections]['group']

        reply = QMessageBox.question(self, 'Confirm Remove', f"Remove item for '{item_text}'?", QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            if row < num_sections:
                self.section_data.pop(row)
            else:
                self.bc_data.pop(row - num_sections)
            self.update_list_view()

    def update_list_view(self):
        """Updates the list view with both section and BC data."""
        display_list = []
        for item in self.section_data:
            params_str = ", ".join([f"{k}={v}" for k, v in item['params'].items()])
            display_list.append(f"[Section] {item['group']}: {item['type']}, {params_str}")
        
        for item in self.bc_data:
            details = ""
            if item['type'] == 'Fix':
                fixes = [axis for axis in ['x', 'y', 'z'] if item.get(f'fix_{axis}')]
                details = f"Fix ({', '.join(fixes).upper() or 'None'})"
            elif item['type'] == 'Force':
                forces = f"F=({item.get('force_x', 0)}, {item.get('force_y', 0)}, {item.get('force_z', 0)})"
                details = f"Force {forces}"
            display_list.append(f"[BC] {item['group']}: {details}")

        self.list_model.setStringList(display_list)

    def mesh_update(self): pass

    def run_simulation(self):
        if not self.mesh:
            QMessageBox.warning(self, "Error", "Please load a mesh file first.")
            return
        
        try:
            E = float(self.young_input.text())
            nu = float(self.poisson_input.text())
            G = E / (2 * (1 + nu))
            pd = 6
            num_nodes = len(self.points)

            props_map = {}
            for sec_data in self.section_data:
                props_map[sec_data['group']] = calculate_section_properties(
                    sec_data['type'], sec_data['params']
                )

            group_id_to_name = {v[0]: k for k, v in self.mesh.field_data.items()}
            line_phys_tags = self.mesh.cell_data_dict['gmsh:physical']['line']

            k_global = np.zeros((pd * num_nodes, pd * num_nodes))
            eps = 1e-6

            for i, element in enumerate(self.conn):
                phys_tag = line_phys_tags[i]
                group_name = group_id_to_name.get(phys_tag)

                if not group_name or group_name not in props_map:
                    QMessageBox.critical(self, "Error", f"Section properties not defined for physical group '{group_name}'.")
                    return
                
                A, I_x, I_y, J, kappa_y, kappa_z, c_y_max, c_z_max = props_map[group_name]

                p1, p2 = self.points[element[0]], self.points[element[1]]
                L = np.linalg.norm(p2 - p1)
                k_ = self.get_timoshenko_stiffness_matrix(L, E, G, A, I_x, I_y, J, kappa_y, kappa_z)
                
                direction = (p2 - p1) / L
                Cxx, Cyx, Czx = direction
                
                if Cxx**2 + Cyx**2 < eps**2:
                    lambda_ = np.array([[0., 0., 1. if Czx > 0 else -1.], [0., 1., 0.], [-1. if Czx > 0 else 1., 0., 0.]])
                else:
                    D = np.sqrt(Cxx**2 + Cyx**2)
                    lambda_ = np.array([[Cxx, Cyx, Czx], [-Cyx/D, Cxx/D, 0], [-Cxx*Czx/D, -Cyx*Czx/D, D]])

                R_theta = np.kron(np.eye(4, dtype=float), lambda_)
                k_local = R_theta.T @ k_ @ R_theta
                
                for j, J_node in enumerate(element):
                    for l, L_node in enumerate(element):
                        k_global[6*J_node:6*J_node+6, 6*L_node:6*L_node+6] += k_local[6*j:6*j+6, 6*l:6*l+6]

            self.u = np.zeros(pd * num_nodes)
            f = np.zeros(pd * num_nodes)
            Up_node = []

            for bc in self.bc_data:
                nodes_to_apply = self.bc_nodes_indexing('vertex', bc['group'])
                for num in nodes_to_apply:
                    if bc['type'] == 'Fix':
                        if bc.get('fix_x'): Up_node.append(6*num + 0)
                        if bc.get('fix_y'): Up_node.append(6*num + 1)
                        if bc.get('fix_z'): Up_node.append(6*num + 2)
                        Up_node.extend([6*num+3, 6*num+4, 6*num+5])
                    elif bc['type'] == 'Force':
                        f[6*num + 0] += bc.get('force_x', 0)
                        f[6*num + 1] += bc.get('force_y', 0)
                        f[6*num + 2] += bc.get('force_z', 0)

            Up_node = sorted(list(set(Up_node)))
            for dof in Up_node: self.u[dof] = 0.0

            Fp_node = [x for x in range(pd*num_nodes) if x not in Up_node]
            k_pu = k_global[np.ix_(Fp_node, Fp_node)]
            k_pp = k_global[np.ix_(Fp_node, Up_node)]
            f_ = f[Fp_node] - k_pp @ self.u[Up_node]
            self.u[Fp_node] = np.linalg.solve(k_pu, f_)
            self.f = k_global @ self.u

            nodal_stresses = np.zeros(num_nodes)
            node_contributions = np.zeros(num_nodes, dtype=int)
            for i, element in enumerate(self.conn):
                phys_tag = line_phys_tags[i]
                group_name = group_id_to_name.get(phys_tag)
                A, I_x, I_y, J, kappa_y, kappa_z, c_y_max, c_z_max = props_map[group_name]

                node1_idx, node2_idx = element
                p1, p2 = self.points[node1_idx], self.points[node2_idx]
                L = np.linalg.norm(p2 - p1)
                
                k_ = self.get_timoshenko_stiffness_matrix(L, E, G, A, I_x, I_y, J, kappa_y, kappa_z)
                
                direction = (p2 - p1) / L
                Cxx, Cyx, Czx = direction
                if Cxx**2 + Cyx**2 < eps**2:
                    lambda_ = np.array([[0.,0.,1. if Czx>0 else -1.],[0.,1.,0.],[-1. if Czx>0 else 1.,0.,0.]])
                else:
                    D = np.sqrt(Cxx**2+Cyx**2)
                    lambda_ = np.array([[Cxx,Cyx,Czx],[-Cyx/D,Cxx/D,0],[-Cxx*Czx/D,-Cyx*Czx/D,D]])
                
                R_theta = np.kron(np.eye(4, dtype=float), lambda_)
                u_global_element = np.concatenate((self.u[6*node1_idx:6*node1_idx+6], self.u[6*node2_idx:6*node2_idx+6]))
                f_local = k_ @ (R_theta @ u_global_element)

                sigma_axial = f_local[6] / A if A > 0 else 0
                sigma_bending1 = abs(f_local[4] * c_z_max / I_x if I_x > 0 else 0) + abs(f_local[5] * c_y_max / I_y if I_y > 0 else 0)
                sigma_bending2 = abs(f_local[10] * c_z_max / I_x if I_x > 0 else 0) + abs(f_local[11] * c_y_max / I_y if I_y > 0 else 0)

                nodal_stresses[node1_idx] += sigma_axial + sigma_bending1
                nodal_stresses[node2_idx] += sigma_axial + sigma_bending2
                node_contributions[node1_idx] += 1
                node_contributions[node2_idx] += 1



            self.smoothed_stresses = np.divide(nodal_stresses, node_contributions, where=node_contributions!=0)
            QMessageBox.information(self, "Success", "Simulation completed.")
  

        except Exception as e:
            QMessageBox.critical(self, "Simulation Error", f"An error occurred during analysis:\n{e}")

    def plot_results(self):
        """Visualizes the analysis results."""
        if self.u is None or self.smoothed_stresses is None:
            QMessageBox.warning(self, "No Data", "Please run the simulation first.")
            return

        pd = 6; scale_factor = 1
        x_orig, y_orig, z_orig = self.points[:,0], self.points[:,1], self.points[:,2]
        x_def = x_orig + scale_factor * self.u[0::pd]
        y_def = y_orig + scale_factor * self.u[1::pd]
        z_def = z_orig + scale_factor * self.u[2::pd]

        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')
        plt.title("3D Timoshenko Beam Analysis")
        ax.set_xlabel("X (m)"); ax.set_ylabel("Y (m)"); ax.set_zlabel("Z (m)")
        

        for n1, n2 in self.conn:
            ax.plot([x_orig[n1], x_orig[n2]], [y_orig[n1], y_orig[n2]], [z_orig[n1], z_orig[n2]], 'k-', lw=1, alpha=0.5, label="Original" if n1==self.conn[0][0] else "")
            ax.plot([x_def[n1], x_def[n2]], [y_def[n1], y_def[n2]], [z_def[n1], z_def[n2]], 'r--', lw=1.5, label="Deformed" if n1==self.conn[0][0] else "")
        

        sc = ax.scatter(x_def, y_def, z_def, c=self.smoothed_stresses, cmap='jet', s=50, edgecolor='k')
        cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.1)
        cbar.set_label('Max Stress (MPa)', rotation=270, labelpad=15)
        clim = sc.get_clim()
        cbar.set_ticks(np.linspace(clim[0], clim[1], 6))
        cbar.set_ticklabels([f'{v/1e6:.1f}' for v in np.linspace(clim[0], clim[1], 6)])
        

        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
        plt.grid(True)

        ax.view_init(elev=20, azim=-60)

        plt.tight_layout()
        plt.show()

    def get_timoshenko_stiffness_matrix(self, L, E, G, A, I_x, I_y, J, kappa_y, kappa_z):
        # phi_z is for bending about z-axis (horizontal), shear force is in y-direction.
        phi_z = (12 * E * I_y) / (G * kappa_y * A * L**2) if (G * kappa_y * A * L**2) > 0 else 0
        # phi_y is for bending about y-axis (vertical), shear force is in z-direction.
        phi_y = (12 * E * I_x) / (G * kappa_z * A * L**2) if (G * kappa_z * A * L**2) > 0 else 0
        
        k11_z = (12*E*I_y)/(L**3*(1+phi_z)) if L > 0 else 0; k12_z = (6*E*I_y)/(L**2*(1+phi_z)) if L > 0 else 0
        k22_z = ((4+phi_z)*E*I_y)/(L*(1+phi_z)) if L > 0 else 0; k23_z = ((2-phi_z)*E*I_y)/(L*(1+phi_z)) if L > 0 else 0
        k11_y = (12*E*I_x)/(L**3*(1+phi_y)) if L > 0 else 0; k12_y = (6*E*I_x)/(L**2*(1+phi_y)) if L > 0 else 0
        k22_y = ((4+phi_y)*E*I_x)/(L*(1+phi_y)) if L > 0 else 0; k23_y = ((2-phi_y)*E*I_x)/(L*(1+phi_y)) if L > 0 else 0
        
        torsion_stiffness = G*J/L if L>0 else 0

        return np.array([
            [A*E/L if L>0 else 0,0,0,0,0,0,-A*E/L if L>0 else 0,0,0,0,0,0],
            [0,k11_z,0,0,0,k12_z,0,-k11_z,0,0,0,k12_z],
            [0,0,k11_y,0,-k12_y,0,0,0,-k11_y,0,-k12_y,0],
            [0,0,0,torsion_stiffness,0,0,0,0,0,-torsion_stiffness,0,0],
            [0,0,-k12_y,0,k22_y,0,0,0,k12_y,0,k23_y,0],
            [0,k12_z,0,0,0,k22_z,0,-k12_z,0,0,0,k23_z],
            [-A*E/L if L>0 else 0,0,0,0,0,0,A*E/L if L>0 else 0,0,0,0,0,0],
            [0,-k11_z,0,0,0,-k12_z,0,k11_z,0,0,0,-k12_z],
            [0,0,-k11_y,0,k12_y,0,0,0,k11_y,0,k12_y,0],
            [0,0,0,-torsion_stiffness,0,0,0,0,0,torsion_stiffness,0,0],
            [0,0,-k12_y,0,k23_y,0,0,0,k12_y,0,k22_y,0],
            [0,k12_z,0,0,0,k23_z,0,-k12_z,0,0,0,k22_z]
        ])

    def bc_nodes_indexing(self, element_type, bc_name):
        """Finds node indices for a given physical group name."""
        try:
            cells = self.mesh.cells_dict.get(element_type)
            if cells is None: return np.array([], dtype=int)
            
            phys_data_map = self.mesh.cell_data_dict.get("gmsh:physical", {})
            phys_data = phys_data_map.get(element_type)
            if phys_data is None: return np.array([], dtype=int)

            target_id = self.mesh.field_data[bc_name][0]
            return np.unique(cells[phys_data == target_id].flatten())
        except (KeyError, IndexError):
            QMessageBox.warning(self, "Warning", f"Physical group '{bc_name}' not found.")
            return np.array([], dtype=int)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = BeamAnalysisWindow()
    window.show()
    sys.exit(app.exec_())