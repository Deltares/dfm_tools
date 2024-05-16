import numpy as np

__all__ = ["LineBuilder"]


class LineBuilder:
    """
    https://matplotlib.org/stable/users/explain/event_handling.html
    """
    def __init__(self, ax):
        # self.waiting_for_entry = True
        line, = ax.plot([], [],'o-') # empty line
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
        self.cid = line.figure.canvas.mpl_connect('key_press_event', self)
        print("Draw a line interactively:\n"
              "- ctrl+leftmouseclick to draw a line\n"
              "- TODO: rightmouseclick to undo\n"
              "- TODO: ctrl+doublemouseclick to wrap up and return line_array\n"
              )
        

    @property
    def line_array(self):
        line_array = np.c_[self.xs, self.ys]
        return line_array
    
    def __call__(self, event):
        # print('click', event)
        
        if not hasattr(event, "button"):
            # don't do anything no mouse button is pressed
            # key event has no button: AttributeError: 'KeyEvent' object has no attribute 'button'
            return
        
        if event.key != "control":
            # don't do anything if ctrl is not pressed
            return
        
        if event.button == 3:
            print('WARNING: right mouse click is ignored')
            return            
        
        if event.dblclick:
            print('WARNING: double click is ignored')
            return
        
        if event.inaxes!=self.line.axes:
            print('WARNING: that click was outside the figure axis, it is ignored')
            return
        
        # otherwise append to line
        print(f"x={event.xdata}m y={event.ydata}")
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        # self.line_array = 
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()

        
# key_press_event: 'canvas', 'guiEvent', 'inaxes', 'key', 'lastevent', 'modifiers', 'name', 'x', 'xdata', 'y', 'ydata']
# button_press_event: 'button', 'canvas', 'dblclick', 'guiEvent', 'inaxes', 'key', 'lastevent', 'modifiers', 'name', 'step', 'x', 'xdata', 'y', 'ydata']