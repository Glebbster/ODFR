import sys
import subprocess

def ensure(pkg, import_name=None):
    try:
        __import__(import_name or pkg)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])

ensure("numpy")
ensure("pandas")
ensure("matplotlib")


import tkinter as tk
from tkinter import filedialog, messagebox
from pathlib import Path

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

read_type = 'fwf' #fwf or csv
SKIPROWS = 56

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("EPR GUI")
        self.geometry("1500x900")
        self.skip_rows = SKIPROWS

        self.default_path = Path(
            "~/Projects/Fortran_prog_for_modeling/data/ODFR4_3_apptest"
        ).expanduser()

        self.dt = None
        self.file_path = None

        left = tk.Frame(self, padx=10, pady=10)
        left.pack(side=tk.LEFT, fill=tk.Y)

        right = tk.Frame(self, padx=10, pady=10, width=100)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        tk.Button(left, text="1) Считать спектр", width=35, command=self.load_spectrum).pack(pady=6)
        tk.Button(left, text="2) Отрисовать", width=35, command=self.plot_simple).pack(pady=6)
        tk.Button(left, text="3) Установить границы (добавить линию)", width=35, command=self.add_border_line).pack(pady=6)
        tk.Button(left,text="4) Отрисовать с границами",width=35,command=self.plot_with_borders_and_coefs).pack(pady=6)


        tk.Label(left, text="Показ кривых:").pack(pady=(12, 4), anchor="w")
        self.show_exp = tk.BooleanVar(value=True)
        self.show_mod = tk.BooleanVar(value=True)
        tk.Checkbutton(left, text="exp", variable=self.show_exp, command=self.update_visibility).pack(anchor="w")
        tk.Checkbutton(left, text="mod", variable=self.show_mod, command=self.update_visibility).pack(anchor="w")

        tk.Label(left, text="Границы: ЛКМ по линии — тащить.\nПКМ по линии — удалить.").pack(
            pady=(12, 0), anchor="w"
        )

        # --- ввод коэффициентов и N ---
        tk.Label(left, text="Коэфициенты по областям слева направо (через запятую):").pack(pady=(14, 4), anchor="w")
        self.coefs_entry = tk.Entry(left, width=38)
        self.coefs_entry.insert(0, "1,2,3")
        self.coefs_entry.pack(anchor="w")

        row = tk.Frame(left)
        row.pack(pady=(8, 0), anchor="w")
        tk.Label(row, text="N точек:").pack(side=tk.LEFT)
        self.npoints_entry = tk.Entry(row, width=8)
        self.npoints_entry.insert(0, "200")
        self.npoints_entry.pack(side=tk.LEFT, padx=(6, 0))

        xlim_row = tk.Frame(left)
        xlim_row.pack(pady=(6, 0), anchor="w")
        tk.Label(xlim_row, text="xlim (xmin,xmax):").pack(side=tk.LEFT)
        self.xlim_entry = tk.Entry(xlim_row, width=16)
        self.xlim_entry.insert(0, "")
        self.xlim_entry.pack(side=tk.LEFT, padx=(6, 0))

        copies_row = tk.Frame(left)
        copies_row.pack(pady=(6, 0), anchor="w")
        tk.Label(copies_row, text="Копий:").pack(side=tk.LEFT)
        self.copies_entry = tk.Entry(copies_row, width=8)
        self.copies_entry.insert(0, "1")
        self.copies_entry.pack(side=tk.LEFT, padx=(6, 0))


        tk.Button(left, text="Сгенерировать weight", width=35, command=self.generate_weights).pack(pady=8)

        self.status = tk.Label(left, text="Файл: (не выбран)\nСтрок: 0\nГраницы: 0", justify="left")
        self.status.pack(pady=12, anchor="w")

        # --- matplotlib ---
        self.fig = Figure(figsize=(12, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar_frame = tk.Frame(right)
        toolbar_frame.pack(fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()

        self.line_exp = None
        self.line_mod = None

        # --- границы (линии) ---
        self.border_lines = []
        self.active_border = None
        self._dragging = False

        self.canvas.mpl_connect("pick_event", self.on_pick)
        self.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.canvas.mpl_connect("button_release_event", self.on_release)
        self.canvas.mpl_connect("button_press_event", self.on_button_press)


    # ---------- чтение ----------
    def load_spectrum(self):
        fp = filedialog.askopenfilename(
            title="Выбери файл спектра",
            initialdir=self.default_path,
            filetypes=[ ("DAT files", "*.dat"),("ESR files", "*.esr"), ("All files", "*.*")]
        )
        if not fp:
            return

        p = Path(fp).expanduser()
        with open(p,'r') as f:
            text = np.asarray(f.readlines())
            # print(text[3])
            idx = 0
            for i,row in enumerate(text):
                if 'Hexp       ,     Sexp     ,     Hmod    ,      Smod' in row:
                    idx = i+1
                    break
        try:
            reader = {
                "csv": lambda: pd.read_csv(
                    p,
                    skiprows=idx,
                    header=None,
                    sep=r"[,\s]+",
                    engine="python",
                    usecols=[0, 1, 2, 3],
                    names=["Hexp", "Sexp", "Hmod", "Smod"],
                    na_values=["", "nan", "NaN"],
                ),
                "fwf": lambda: pd.read_fwf(
                    p,
                    colspecs=[(0,15),(16,34),(35,50),(51,66),(67,88),(89,107)],
                    skiprows=idx,
                    header=None,
                    usecols=[0,1,2,3,4,5],
                    names=["Hexp","Sexp","Hmod","Smod","Smod2","Smod3"],
                    na_values=[",","", "nan", "NaN"],
                ),
            }[read_type]

            dt = reader()

        except Exception:
            try:
                raw = pd.read_csv(
                    p,
                    skiprows=idx,
                    header=None,
                    sep=r"[,\s]+",
                    engine="python",
                    na_values=["", "nan", "NaN"],
                )

                raw = raw.dropna(how="all")

                while raw.shape[1] < 4:
                    raw[raw.shape[1]] = np.nan

                dt = raw.iloc[:, :4]
                dt.columns = ["Hexp", "Sexp", "Hmod", "Smod"]

            except Exception as e:
                messagebox.showerror("Ошибка чтения", repr(e))
                return

        self.dt = dt
        self.file_path = p

        # сброс границ при загрузке нового файла
        self._clear_borders()

        self.status.config(text=f"Строк: {len(dt)}\nГраницы: {len(self.border_lines)}")

    # ---------- plot ----------
    def plot_simple(self):
        if self.dt is None:
            messagebox.showwarning("Нет данных", "Сначала нажми: 1) Считать спектр")
            return

        self.ax.clear()
        self.ax.grid(True)

        self.line_exp, = self.ax.plot(self.dt["Hexp"], self.dt["Sexp"], label="exp")
        self.line_mod, = self.ax.plot(self.dt["Hmod"], self.dt["Smod"], label="mod")
        txt = self.xlim_entry.get().strip()
        if txt:
            try:
                xmin, xmax = [float(v) for v in txt.split(",")]
                self.ax.set_xlim(xmin, xmax)
            except Exception:
                pass

        # перерисовать границы поверх графика (если были)
        old_xs = sorted([float(ln.get_xdata()[0]) for ln in self.border_lines])
        self._clear_borders(redraw=False)
        for x in old_xs:
            self.border_lines.append(self._make_border_line(x))

        self.update_visibility()
        self.ax.legend()
        self.canvas.draw()
        self.status.config(text=f"Строк: {len(self.dt)}\nГраницы: {len(self.border_lines)}")
        self.ax.tick_params(labelsize = 28)

    def update_visibility(self):
        if self.line_exp is not None:
            self.line_exp.set_visible(self.show_exp.get())
        if self.line_mod is not None:
            self.line_mod.set_visible(self.show_mod.get())
        self.canvas.draw_idle()

    # ---------- границы ----------
    def _make_border_line(self, x: float):
        ln = self.ax.axvline(x, color="red", linewidth=2, alpha=0.9)
        ln.set_picker(7)
        return ln

    def _highlight_active(self):
        for ln in self.border_lines:
            ln.set_linewidth(3 if ln is self.active_border else 2)
            ln.set_alpha(1.0 if ln is self.active_border else 0.7)

    def _clear_borders(self, redraw=True):
        for ln in self.border_lines:
            try:
                ln.remove()
            except Exception:
                pass
        self.border_lines.clear()
        self.active_border = None
        if redraw:
            self.canvas.draw_idle()

    def add_border_line(self):
        if self.dt is None:
            messagebox.showwarning("Нет данных", "Сначала нажми: 1) Считать спектр")
            return

        xmin, xmax = self.ax.get_xlim()
        x = 0.5 * (xmin + xmax)
        ln = self._make_border_line(x)
        self.border_lines.append(ln)

        self.active_border = ln
        self._highlight_active()
        self.canvas.draw_idle()
        self.status.config(text=f"Строк: {len(self.dt)}\nГраницы: {len(self.border_lines)}")

    def on_pick(self, event):
        artist = event.artist
        if artist in self.border_lines:
            self.active_border = artist
            self._highlight_active()
            self._dragging = True
            self.canvas.draw_idle()

    def on_motion(self, event):
        if not self._dragging or self.active_border is None:
            return
        if event.inaxes != self.ax or event.xdata is None:
            return
        x = float(event.xdata)
        self.active_border.set_xdata([x, x])
        self.canvas.draw_idle()

    def on_release(self, _event):
        self._dragging = False

    def on_button_press(self, event):
        if event.inaxes != self.ax:
            return
        if event.button != 3:
            return

        target = None
        for ln in self.border_lines:
            hit, _ = ln.contains(event)
            if hit:
                target = ln
                break

        if target is not None:
            if target is self.active_border:
                self.active_border = None
            target.remove()
            self.border_lines.remove(target)
            self._highlight_active()
            self.canvas.draw_idle()
            if self.dt is not None:
                self.status.config(text=f"Строк: {len(self.dt)}\nГраницы: {len(self.border_lines)}")


    # ---------- weights ----------
    def generate_weights(self):
        if self.dt is None:
            messagebox.showwarning("Нет данных", "Сначала нажми: 1) Считать спектр")
            return

        # N
        try:
            npoints = int(self.npoints_entry.get().strip())
            if npoints <= 1:
                raise ValueError
        except Exception:
            messagebox.showerror("Ошибка", "N точек должно быть целым числом > 1")
            return
        # copies
        try:
            copies = int(self.copies_entry.get().strip())
            if copies < 1:
                raise ValueError
        except Exception:
            messagebox.showerror("Ошибка", "Копий должно быть целым числом >= 1")
            return
        

        # коэффициенты
        try:
            coefs = [float(x.strip()) for x in self.coefs_entry.get().split(",") if x.strip() != ""]
        except Exception:
            messagebox.showerror("Ошибка", "Коэфициенты: только числа через запятую")
            return

        # границы по X: xmin + (красные линии) ; области = len(lines)+1
        xmin = float(self.dt["Hexp"].min())
        xmax = float(self.dt["Hexp"].max())

        lines_x = sorted([float(ln.get_xdata()[0]) for ln in self.border_lines])

        starts = np.array([xmin] + lines_x, dtype=float)

        if len(coefs) != len(starts):
            messagebox.showerror(
                "Ошибка",
                f"Нужно коэфициентов: {len(starts)} (областей = границы+1). Сейчас: {len(coefs)}.\n"
                f"Границ: {len(lines_x)}"
            )
            return

        coefs_arr = np.array(coefs, dtype=float)

        numspace = np.linspace(xmin, xmax, npoints)

        idx = np.searchsorted(starts, numspace, side="right") - 1
        idx = np.clip(idx, 0, len(coefs_arr) - 1)
        weights = coefs_arr[idx]
        if copies > 1:
            weights = np.tile(weights, copies)
        

        try:
            np.savetxt("weight", weights.reshape(-1, 1), fmt="%.1f")
        except Exception as e:
            messagebox.showerror("Ошибка записи", repr(e))
            return

        messagebox.showinfo(
            "Ок",
            f"weights сохранен\nN={len(weights)}\nstarts={starts.tolist()}\ncoefs={coefs}"
        )
    def plot_with_borders_and_coefs(self):
        if self.dt is None:
            messagebox.showwarning("Нет данных", "Сначала считай спектр")
            return

        # парсинг коэффициентов
        try:
            coefs = [float(x.strip()) for x in self.coefs_entry.get().split(",") if x.strip()]
        except Exception:
            messagebox.showerror("Ошибка", "Коэфициенты должны быть числами через запятую")
            return

        # границы
        xmin = float(self.dt["Hexp"].min())
        xmax = float(self.dt["Hexp"].max())
        lines_x = sorted([float(ln.get_xdata()[0]) for ln in self.border_lines])
        starts = [xmin] + lines_x

        if len(coefs) != len(starts):
            messagebox.showerror(
                "Ошибка",
                f"Коэфициентов должно быть {len(starts)} (областей)"
            )
            return

        # базовый график
        self.plot_simple()

        # y для подписей
        y_text = self.ax.get_ylim()[1]


        # подписи коэффициентов
        for i, coef in enumerate(coefs):
            x0 = starts[i]
            x1 = starts[i + 1] if i < len(starts) - 1 else xmax
            x_mid = 0.5 * (x0 + x1)

            self.ax.annotate(
            str(coef),
            xy=(x_mid, y_text),
            xytext=(0, 8),               # ↑ пиксели НАД графиком
            textcoords="offset points",
            ha="center",
            va="bottom",
            color="red",
            fontsize=28,
            fontweight="bold",
            clip_on=False                # ← важно
            )


        self.canvas.draw_idle()



if __name__ == "__main__":
    App().mainloop()
