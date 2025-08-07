import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from scipy import stats
import frc
import math
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class FRCAnalyzer:
    def __init__(self, master):
        self.master = master
        self.master.title("SMLM FRC Analyzer")
        self.master.protocol("WM_DELETE_WINDOW", self.master.quit)

        # GUI layout frame
        top_frame = tk.Frame(master)
        top_frame.pack(pady=5)

        # File selection
        self.browse_button = tk.Button(top_frame, text="Browse CSV", command=self.browse_file)
        self.browse_button.pack(padx=5)

        self.file_label = tk.Label(top_frame, text="No file selected", width=70, justify="left")
        self.file_label.pack(anchor="w", padx=5)

        # Pixel size entry
        pixel_frame = tk.Frame(master)
        pixel_frame.pack(pady=5)

        self.pixel_size_label = tk.Label(pixel_frame, text="Pixel size (nm):")
        self.pixel_size_label.pack(side=tk.LEFT, padx=5)

        self.pixel_size_entry = tk.Entry(pixel_frame, width=10)
        self.pixel_size_entry.insert(0, "10")
        self.pixel_size_entry.pack(side=tk.LEFT)
        self.pixel_size_entry.bind("<Return>", lambda event: self.update_images())

        self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(10, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().pack()

        self.selector = None
        self.spots = None
        self.img1 = None
        self.img2 = None
        self.last_roi = None

    def browse_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if file_path:
            self.file_label.config(text=file_path)
            self.load_data(file_path)

    def load_data(self, path):
        self.spots = np.genfromtxt(path, delimiter=",", names=True)
        self.update_images()

    def update_images(self):
        try:
            pixel_size = float(self.pixel_size_entry.get())
        except ValueError:
            print("Invalid pixel size")
            return

        half1 = self.spots[self.spots['frame'] % 2 == 0]
        half2 = self.spots[self.spots['frame'] % 2 == 1]

        image_size_auto = np.max([np.max(half1['x_nm']), np.max(half2['x_nm']), np.max(half1['y_nm']), np.max(half2['y_nm'])])
        image_size_auto_rounded = math.ceil((1.1*image_size_auto)/100) * 100
        image_size_nm = (image_size_auto_rounded, image_size_auto_rounded)
        self.image_size_px = tuple(int(s / pixel_size) for s in image_size_nm)

        self.img1 = self.render_image(half1, pixel_size)
        self.img2 = self.render_image(half2, pixel_size)

        self.ax1.clear()
        self.ax1.imshow(self.img1, cmap='hot', vmax=np.percentile(self.img1, 99))
        self.ax1.set_title("Half 1")

        if self.selector:
            self.selector.set_active(False)

        self.selector = RectangleSelector(self.ax1, self.onselect, useblit=True, button=[1], interactive=True, spancoords='pixels')
        self.ax2.clear()
        self.ax2.set_title("FRC")
        self.ax2.set_xlabel("Spatial frequency (1/nm)")
        self.canvas.draw()

        if self.last_roi:
            self.calculate_frc(*self.last_roi)

    def render_image(self, spots_subset, pixel_size):
        image = np.zeros(self.image_size_px)
        for spot in spots_subset:
            x_px = int(spot['x_nm'] / pixel_size)
            y_px = int(spot['y_nm'] / pixel_size)
            intensity = int(spot['intensity_photon'])
            if 0 <= x_px < self.image_size_px[0] and 0 <= y_px < self.image_size_px[1]:
                image[y_px, x_px] += 1
        sigma_nm = stats.mode(spots_subset['uncertainty_xy_nm'], keepdims=True).mode[0]
        sigma_px = sigma_nm / pixel_size
        return gaussian_filter(image, sigma=sigma_px)

    def onselect(self, eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        width = abs(x2 - x1)
        height = abs(y2 - y1)
        side = min(width, height)
        x0 = min(x1, x2)
        y0 = min(y1, y2)

        self.last_roi = (int(x0), int(y0), int(side), int(side))
        self.calculate_frc(*self.last_roi)

    def calculate_frc(self, x, y, w, h):
        pixel_size = float(self.pixel_size_entry.get())

        x1 = math.ceil(x/10)*10
        y1 = math.ceil(y/10)*10
        w1 = math.floor(w/10)*10
        h1 = math.floor(h/10)*10

        img1_cropped = self.img1[y1:y1+h1, x1:x1+w1]
        img2_cropped = self.img2[y1:y1+h1, x1:x1+w1]

        img1_cropped = frc.util.square_image(img1_cropped, add_padding=False)
        img2_cropped = frc.util.square_image(img2_cropped, add_padding=False)
        img1_cropped = frc.util.apply_tukey(img1_cropped)
        img2_cropped = frc.util.apply_tukey(img2_cropped)

        frc_curve = frc.two_frc(img1_cropped, img2_cropped)
        smoothed_frc = savgol_filter(frc_curve, window_length=15, polyorder=2)

        img_size = img1_cropped.shape[0]
        xs_pix = np.arange(len(frc_curve)) / img_size
        scale = 1 / pixel_size
        xs_nm_freq = xs_pix * scale

        try:
            frc_res, res_y, thres = frc.frc_res(xs_nm_freq, smoothed_frc, img_size)
        except frc.deps_types.NoIntersectionException:
            messagebox.showwarning("FRC Error", "No intersection found in FRC curve.")
            return

        self.ax2.clear()
        self.ax2.plot(xs_nm_freq, thres(xs_nm_freq), label='Threshold')
        self.ax2.plot(xs_nm_freq, smoothed_frc, label='FRC')
        self.ax2.axvline(x=1/frc_res, color='green', linestyle='--', label=f"Resolution = {frc_res:.1f} nm")
        self.ax2.set_title("FRC")
        self.ax2.set_xlabel("Spatial frequency (1/nm)")
        self.ax2.legend()
        self.canvas.draw()

        print(f"Estimated resolution: {frc_res:.1f} nm")

if __name__ == '__main__':
    root = tk.Tk()
    app = FRCAnalyzer(root)
    root.mainloop()
