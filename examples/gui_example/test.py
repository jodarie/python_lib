from Tkinter import *
from datetime import datetime
import math

BOX_WIDTH = 330
BOX_HEIGHT = 200

class Application(Frame):

    def ap_canvas(self):
        self.canvas = Canvas(root, width = BOX_WIDTH, height = BOX_HEIGHT, bg = 'white')
        self.canvas.pack()
        
    def quit_button(self):
        self.QUIT = Button(self)
        self.QUIT["text"] = "QUIT"
        self.QUIT["fg"] = "red"
        self.QUIT["command"] = self.quit
        self.QUIT.pack({"side": "left"})
        
    def time_button(self):
        self.time = Button(self)
        self.time["text"] = "Time"
        self.time["command"] = self.print_time
        self.time.pack({"side": "left"})

    def sq_root_button(self):
        self.sq_root = Button(self)
        self.sq_root ["text"] = "Square Root"
        self.sq_root ["command"] = self.print_sq_root
        self.sq_root.pack({"side": "left"})

    def entry_box(self):
        self.box = Entry(self)
        self.box.pack({"side": "left"})
        self.box.bind('<Return>', self.get_input)

    def print_time(self):
        text = 'The current time is '+str(datetime.now().time())
        self.print_to_canvas(text)

    def print_sq_root(self):
        if self.input is not None:
            text = 'The square root of ' + str(self.input) + ' is ' + str(math.sqrt(self.input)) + '.'
        else:
            text = 'Error: Input value not specified.'
        self.print_to_canvas(text)

    def print_to_canvas(self, text):
        self.canvas.delete('all')
        self.canvas.create_text(10, 50 , font = 'Verdana 12', text = text, anchor = 'w')

    def get_input(self, input):
        self.input = float(self.box.get())
        self.box.delete(0,END)
        text = 'The input value is ' + str(self.input) + '.'
        self.print_to_canvas(text)

    def __init__(self, master=None):
        self.input = None
        Frame.__init__(self, master)
        self.pack()
        self.quit_button()
        self.time_button()
        self.sq_root_button()
        self.entry_box()
        self.ap_canvas()

root = Tk()
root.title('GUI TEST')
app = Application(master=root)
app.mainloop()
root.destroy()
