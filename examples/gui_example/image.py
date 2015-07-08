from Tkinter import *

root = Tk()

canvas = Canvas(root)

canvas.grid(row = 0, column = 0)

photo = PhotoImage(file = './flags/botswana.gif')

canvas.create_image(0,0, image=photo)
    
root.mainloop()
