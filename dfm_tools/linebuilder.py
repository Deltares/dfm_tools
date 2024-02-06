import numpy as np

__all__ = ["LineBuilder"]


class LineBuilder:
    """
    https://matplotlib.org/stable/users/explain/event_handling.html
    """
    def __init__(self, line):
        #self.waiting_for_entry = True
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        print('click', event)
        # if event.inaxes!=self.line.axes:
        #     print('WARNING: that click was outside the figure axis')
        #     return
        # if event.dblclick:
        #     print('WARNING: double click is interpreted as single click')
        #     return
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line_array = np.c_[self.xs, self.ys]
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()
