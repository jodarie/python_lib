#$Id: buttonbox.py,v 1.2 2004/03/17 04:29:31 mandava Exp $
# this is a program to create a button box which contains more than one
# button.
                                                                                                                               
title='Pmw.button box'
                                                                                                                               
#import Pmw from this directory tree
import sys
sys.path[:0] = ['../../..']
                                                                                                                               
import Tkinter
import Pmw
                                                                                                                               
                                                                                                                               
# ButtonBox inherits from pmw.MegaWidget
# labelpos specifies where to place the label component.if the label pos
# is not none it means that a 'frame' component has been created to act as
# a container of the buttons created.the label component should be a
# concatination of one or two of the letters N,S,W,E. in this case the
# label component is 'nw'which means the label would be placed on the top
# of the left hand side.
# if there is no label component ,then no frame component is created and
# the 'hull' component acts as the container.the hull acts as the body
#for the entire megawidget.other components act as children of hull to
#further specialize the widget.
#the label_text option is used to set the text option of the label.
                                                                                                                               
                                                                                                                               
# this class creates a manager widget for containing buttons.                                                                                                                                
class Demo:
     def __init__(self,parent):
  # create and pack the buttonbox
		self.buttonbox = Pmw.ButtonBox(parent,
                labelpos='nw',
                label_text='ButtonBox:')
                                                                                                                               
 #padx and pady specify the padding distance to leave between each button
 #and also between the buttons and the outer edge of the button box.
		self.buttonbox.pack(fill='both', padx=8, pady=8)
 #add some button to the buttonbox
 #'add' method adds a button to the end of the buttonbox.
 # it creates a button with the text keyword argument, in this
 # case three buttons with text Ok,cancel, and apply are created.
                                                                                                                               
		self.buttonbox.add('OK')
		self.buttonbox.add('apply')
		self.buttonbox.add('cancel')
######################################################################
                                                                                                                               
#create demo in root window for testing
                                                                                                                               
if __name__=='__main__':
    root= Tkinter.Tk()
    Pmw.initialise(root)
    root.title(title)
                                                                                                                               
    widget= Demo(root)
    root.mainloop()
                                                                                                                               
                                                                                                                               

