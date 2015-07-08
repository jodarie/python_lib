#$Id: entryfield.py,v 1.2 2004/03/17 04:29:31 mandava Exp $
#this is a program depicting the use of entry field widget.
#entry widgets are basic widgets used to collect input from the user.
#entry widgets are limited to a single line of text which can be in only
#one font. 
#the root is also packed with 4 buttons along with the entry widget..

import Tkinter 
from Tkinter import *
root =Tk()
root.title('entry widget')
Label (text='enter the file name').pack(side=TOP,padx=10,pady=10)
Entry(root, width=10).pack(side=TOP,padx=10,pady=10)
Button(root, text='open').pack(side= LEFT)
Button(root, text='edit').pack(side= LEFT)
Button(root, text='exit').pack(side= RIGHT)
Button(root, text='close').pack(side= RIGHT)
root.mainloop()

