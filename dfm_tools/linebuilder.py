import numpy as np
import matplotlib.pyplot as plt

__all__ = ["LineBuilder"]


class LineBuilder:
    """
    To interactively draw a line in a figure axis, for instance to use as cross-section.
    """
    def __init__(self, ax=None):
        # print some user help
        print("Draw a line interactively:\n"
              "- ctrl+leftmouseclick to add a point to the line\n"
              "- ctrl+rightmouseclick to remove the last point of the line\n"
              "- ctrl+doublemouseclick to finish and let the script continue\n"
              )
        
        # get current axis if not provided
        if ax is None:
            ax = plt.gca()
        
        # initialize x/y arrays and line object
        self.xs = []
        self.ys = []
        line, = ax.plot(self.xs, self.ys,'o-') # empty line
        self.line = line
        
        # register both button press events and key press events
        self.cid_button = self.line.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_key = self.line.figure.canvas.mpl_connect('key_press_event', self.on_press)
        # start a blocking event loop to prevent continuation of the script where LineBuilder was called
        self.line.figure.canvas.start_event_loop()
    
    @property
    def line_array(self):
        line_array = np.c_[self.xs, self.ys]
        return line_array
    
    def add_xy_to_line(self, event):
        print(f"adding point: x={event.xdata:.6f}, y={event.ydata:.6f}")
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()

    def remove_last_xy_from_line(self, event):
        print("removing last point if present")
        self.xs = self.xs[:-1]
        self.ys = self.ys[:-1]
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()

    def finish_linebuilder(self):
        # disconnect the interactive line drawing in the figure
        self.line.figure.canvas.mpl_disconnect(self.cid_button)
        self.line.figure.canvas.mpl_disconnect(self.cid_key)
        # let the rest of the script where LineBuilder was called continue
        self.line.figure.canvas.stop_event_loop()
        print("interactive line drawing finished")
    
    def on_press(self, event):
        # do nothing if ctrl is not pressed
        if event.key != "control":
            return
        
        # do nothing if no mouse button is pressed
        if not hasattr(event, "button"):
            # key event has no button: AttributeError: 'KeyEvent' object has no attribute 'button'
            return
        
        # do nothing if mouse click is outside axis
        if event.inaxes!=self.line.axes:
            print('clicks outside of the figure axis are ignored')
            return
        
        # add new point to line and wrap up line upon "ctrl + left mouse double click"
        if event.dblclick:
            self.finish_linebuilder()
            return
        
        # add new point to line upon "ctrl + left mouse click"
        if event.button == 1:
            self.add_xy_to_line(event)
        
        # remove last point from line upon "ctrl + right mouse click"
        if event.button == 3:
            self.remove_last_xy_from_line(event)
            return            
