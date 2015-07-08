from Tkinter import *
import math, random, numpy as np

BOX_WIDTH = 320
BOX_HEIGHT = 200

IMAGE_DIR = "./flags/"
IMAGE_FILE = "image_list.txt"

FONT = 'Verdana 12'

class Application(Frame):

    def ap_canvas(self):
        self.canvas = Canvas(root, width = BOX_WIDTH,
                             height = BOX_HEIGHT, bg = 'white')
        self.canvas.pack()
        
    def quit_button(self):
        self.QUIT = Button(self)
        self.QUIT["text"] = "QUIT"
        self.QUIT["fg"] = "red"
        self.QUIT["command"] = self.quit
        self.QUIT.pack({"side": "left"})

    def next_button(self):
        self.NEXT = Button(self)
        self.NEXT["text"] = "NEXT"
        self.NEXT["fg"] = "blue"
        self.NEXT["command"] = self.goto_next
        self.NEXT.pack({"side": "left"})

    def new_button(self):
        self.NEW = Button(self)
        self.NEW["text"] = "NEW GAME"
        self.NEW["fg"] = "blue"
        self.NEW["command"] = self.new_game
        self.NEW.pack({"side": "left"})

    def entry_box(self):
        self.box = Entry(self)
        self.box.pack({"side": "left"})
        self.box.bind('<Return>', self.get_input)

    def get_input(self, input):
        input = self.box.get()
        self.box.delete(0,END)
        self.check_answer(input)

    def read_image_file(self):
        self.flag_list = np.genfromtxt(IMAGE_DIR + IMAGE_FILE,
                                       dtype = 'S', unpack = True, delimiter = ',')
        self.flag_number = np.array(self.flag_list[0, :], dtype = 'int') - 1
        self.n_flags = self.flag_list.shape[1]
        
    def show_image(self, i):
        self.canvas.delete('all')
        self.photo = PhotoImage(file = IMAGE_DIR + self.flag_list[2, i])
        self.canvas.create_image(25, 10, image = self.photo, anchor = NW)

    def new_game(self):
        self.game_status = 'active'
        self.score = 0
        self.lives = 3
        random.shuffle(self.flag_number)
        self.goto_next()

    def goto_next(self):
        if self.game_status == 'active':
            if self.count == None:
                self.count = 0
            if self.count == len(self.flag_number) - 1:
                self.count = 0
            else:
                self.count += 1
            self.current_flag = self.flag_number[self.count]
            self.show_image(self.current_flag)
            self.print_score()

    def check_answer(self, answer):
        if self.game_status == 'active':
            if answer == self.flag_list[1, self.current_flag]:
                self.score += 1
                if len(self.flag_number) == 1:
                    self.final_score()
                else:
                    self.flag_number = np.delete(self.flag_number, self.count)
                    self.count = 0
                    self.goto_next()                        
            else:
                self.lives -= 1
                text = 'Incorrect! '
                self.print_to_canvas(text)
                self.print_score()
                if self.lives == 0:
                    self.final_score()

    def print_to_canvas(self, text):
        self.canvas.delete(self.current_guess)
        self.current_guess = self.canvas.create_text(25, 170 , font = FONT, text = text, anchor = 'w')

    def print_score(self):
        self.canvas.delete(self.current_lives)
        text1 = str(self.score) + ' of ' + str(self.n_flags)
        text2 = 'Lives: ' + str(self.lives)
        self.canvas.create_text(235, 50 , font = FONT, text = text1, anchor = 'w')
        self.current_lives = self.canvas.create_text(235, 70 , font = FONT, text = text2, anchor = 'w')

    def final_score(self):
        self.game_status = 'over'
        self.canvas.delete('all')
        text1 = 'You identified ' + str(self.score / self.n_flags * 100) + '% of the flags.'
        text2 = 'Lives Remainig: ' + str(self.lives)
        text3 = 'GAME OVER!'
        self.canvas.create_text(10, 25 , font = FONT, text = text1, anchor = 'w')
        self.canvas.create_text(10, 45 , font = FONT, text = text2, anchor = 'w')
        self.canvas.create_text(110, 95 , font = FONT, text = text3, anchor = 'w')
        
    def __init__(self, master=None):
        self.current_lives = None
        self.flag_list = None
        self.flag_number = None
        self.n_flags = None
        self.count = 0
        self.current_flag = None
        self.current_guess = None
        Frame.__init__(self, master)
        self.pack()
        self.quit_button()
        self.next_button()
        self.new_button()
        self.entry_box()
        self.ap_canvas()
        self.read_image_file()
        self.new_game()

root = Tk()
root.title('WORLD FLAGS')
app = Application(master=root)
app.mainloop()
root.destroy()
