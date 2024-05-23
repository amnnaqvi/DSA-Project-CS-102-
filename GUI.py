from tkinter import *
from SearchADT import*

bg_image = None #bg image is a global variable
root = None #main root
Choose_existing = None #open_Choose_Existing root
addOwn = None #open_Add_Your_Own root
displayChain = None #display_PolyPeptide root
displayStrand = None #open_Display_List root
editStrand = None #open Input_Pattern root
searchroot = None #search root

mutation_root = None #Mutate root
compare_root = None #compare root
sequence, bwt, suffix_array, c_table, occ_ranks, coordinates = None, None, None, None, None, [] #variables for FM index
mutation_count = 0 #allow 3 mutations in a protein before ending by default
first_sequence = None

def fill_dictionary(AA): #fixed colours for the amino acids; TC is O(1)
    color_dict = {'Ala': 'red','Arg': 'green','Asn': 'blue','Asp': 'pink','Cys': 'yellow','Gln': 'orange','Glu': 'purple','Gly': 'cyan','His': 'magenta','Ile': 'bisque','Leu': 'brown','Lys': 'gray','Met': 'lime','Phe': 'turquoise','Pro': 'olive','Ser': 'teal','Thr': 'indigo','Trp': 'violet','Tyr': 'salmon','Val': 'lavender'}
    return color_dict[AA]

def draw_circle_with_text(canvas, x, y, radius, text): #create circles to display polypeptide chain; TC is O(1)
    canvas.create_oval(x - radius, y - radius, x + radius, y + radius, outline="black", fill=fill_dictionary(text))
    canvas.create_text(x, y, text=text)

def main(): #title screen; O(1)

    # Create the root window
    global root, bg_image
    root = Tk()
    root.geometry("1920x1080")
    root.title("DSA Team 7 - Amn & Rayyan's DNA Sequencer")

    # Check if bg_image is None, if so, load the image
    if bg_image is None:
        bg_image = PhotoImage(file="background_image.png")

    # Use a label to display the background image
    bg_label = Label(root, image=bg_image)
    bg_label.place(relx=0, rely=0, relwidth=1, relheight=1)

    # Create the title label for the main page
    title_label = Label(root, text="Gene Searching and Mutations Via FM Index!", font=("Times New Roman", 20, "bold", "underline"), bg="purple", fg="white", borderwidth=3, relief=FLAT)
    title_label.place(relx=0.5, rely=0.2, anchor=CENTER)

    # Create the buttons for the main page
    button1 = Button(root, text="Choose Protein From Library", font=("Times New Roman", 16), width=24, height=3, command=lambda: open_Choose_Existing(), bg="bisque", fg="black", borderwidth=2, relief=RAISED)
    button1.place(relx=0.3, rely=0.55, anchor=CENTER)

    button2 = Button(root, text="Add Your Own", font=("Times New Roman", 16), width=24, height=3, command=lambda: open_Add_Your_Own(), bg="bisque", fg="black", borderwidth=2, relief=RAISED)
    button2.place(relx=0.7, rely=0.55, anchor=CENTER)

    disclaimer_label = Label(root, text="Note: For simplicity, the program will only consider max 243 bases!", font=("Times New Roman", 14, "bold"), bg="purple", fg="white", relief=FLAT)
    disclaimer_label.place(relx=0.5, rely=0.9, anchor=CENTER)

    root.mainloop()

def open_Choose_Existing(): #gives user 10 proteins to choose from; O(1)

    global root, bg_image, Choose_existing

    Choose_existing = Toplevel()
    Choose_existing.geometry("1920x1080")
    Choose_existing.title("Choosing Existing Protein")

    bg_label = Label(Choose_existing, image=bg_image)
    bg_label.place(relx=0, rely=0, relwidth=1, relheight=1)

    title_label = Label(Choose_existing, text="Choose Any One Protein From Existing Library:", font=("Times New Roman", 20, "underline"), bg="black", fg="white", borderwidth=3, relief=SUNKEN)
    title_label.place(anchor=CENTER, relx=0.5, rely=0.1)

    # Create buttons for all proteins
    button_A = Button(Choose_existing, text="Hemoglobin", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Hemoglobin",1), bg="red", fg="white", relief=RAISED)
    button_B = Button(Choose_existing, text="Insulin", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Insulin",1), bg="bisque", fg="black" ,relief=RAISED)
    button_C = Button(Choose_existing, text="Collagen", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Collagen",1), bg="green", fg="white",relief=RAISED)
    button_D = Button(Choose_existing, text="Actin", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Actin",1), bg="blue", fg="white",relief=RAISED)
    button_E = Button(Choose_existing, text="Myosin", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Myosin",1), bg="orange", fg="black",relief=RAISED)
    button_F = Button(Choose_existing, text="Fibrinogen", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Fibrinogen",1), bg="yellow", fg="black",relief=RAISED)
    button_G = Button(Choose_existing, text="Melanin", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Melanin",1), bg="purple", fg="white",relief=RAISED)
    button_H = Button(Choose_existing, text="Glucagon", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Glucagon",1), bg="cyan", fg="black",relief=RAISED)
    button_I = Button(Choose_existing, text="Oxytocin", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Oxytocin",1), bg="pink", fg="black",relief=RAISED)
    button_J = Button(Choose_existing, text="Immunoglobulin G", font=("Times New Roman", 16), width=15, height=2, command=lambda: display_PolyPeptide("Immunoglobulin G",1), bg="brown", fg="white",relief=RAISED)

    button_A.place(relx=0.3, rely=0.25, anchor=CENTER)
    button_B.place(relx=0.5, rely=0.25, anchor=CENTER)
    button_C.place(relx=0.7, rely=0.25, anchor=CENTER)
    button_D.place(relx=0.3, rely=0.45, anchor=CENTER)
    button_E.place(relx=0.5, rely=0.45, anchor=CENTER)
    button_F.place(relx=0.7, rely=0.45, anchor=CENTER)
    button_G.place(relx=0.3, rely=0.65, anchor=CENTER)
    button_H.place(relx=0.5, rely=0.65, anchor=CENTER)
    button_I.place(relx=0.7, rely=0.65, anchor=CENTER)
    button_J.place(relx=0.5, rely=0.85, anchor=CENTER)

    root.withdraw()

    #create a button to return to previous page

    def back():
        Choose_existing.destroy()
        root.deiconify()

    BackButton = Button(Choose_existing, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back, bg="red", fg="white", relief=RAISED)
    BackButton.place(relx=0.1, rely=0.1, anchor=CENTER)

def open_Add_Your_Own(): #allows user to enter Template or Coding strand; O(n)

    global root, bg_image, addOwn

    addOwn = Toplevel()
    addOwn.geometry("1920x1080")
    addOwn.title("Add Your Own Protein")

    bg_label = Label(addOwn, image=bg_image)
    bg_label.place(relx=0, rely=0, relwidth=1, relheight=1)

    root.withdraw()
            
    def process_AddOwn(Protein, num):
        Protein = Protein.upper()
        if num == 2: #convert all coding strands to template for standardization
            Protein = coding_to_template(Protein) #O(n)
        result = validation_check(Protein) #O(n)
        if result == "Processing...":
            choose.config(text="Please Enter Template Strand or Coding Strand:", bg='bisque', fg='black') #reset the title incase the user comes back to this page
            display_PolyPeptide(Protein, 2)
        else:
            choose.config(text=result, fg='red')

    choose = Label(addOwn, text="Please Enter Template Strand or Coding Strand:", font=("Times New Roman", 18, "bold"), width=50, height=1, bg='bisque', fg='black', relief=FLAT)
    choose.place(relx=0.5, rely=0.4, anchor=CENTER)

    OwnSeq = Entry(addOwn, width=60, relief=SUNKEN)
    OwnSeq.place(relx=0.5, rely=0.5, anchor="center")

    validate1_button = Button(addOwn, text="Validate My Template Strand", command=lambda: process_AddOwn(OwnSeq.get(), 1))
    validate1_button.place(relx=0.4, rely=0.6, anchor="center")

    validate2_button = Button(addOwn, text="Validate My Coding Strand", command=lambda: process_AddOwn(OwnSeq.get(), 2))
    validate2_button.place(relx=0.6, rely=0.6, anchor="center")

    def back():
        addOwn.destroy()
        root.deiconify()

    BackButton = Button(addOwn, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back, bg="red", fg="white", relief=RAISED)
    BackButton.place(relx=0.1, rely=0.1, anchor=CENTER)

def display_PolyPeptide(Protein, Option): #display circles representing amino acids; O(n**2)

    global displayChain, Choose_existing, bg_image, addOwn

    displayChain = Toplevel()
    displayChain.geometry("1920x1080")
    displayChain.title("Current Amino Acid Chain")

    def back_1(): #go back to choosing from 10 options
        displayChain.destroy()
        Choose_existing.deiconify()

    def back_2(): #go back to entering own
        displayChain.destroy()
        addOwn.deiconify()

    bg_label = Label(displayChain, image=bg_image)
    bg_label.place(relx=0, rely=0)

    canvas = Canvas(displayChain, width=680, height=680) #creates a white canvas to draw circles on
    canvas.place(anchor=CENTER, relx=0.5,rely=0.5)
    
    if Option == 1:
        gene = csv_reader(Protein) #O(n+m)
        Choose_existing.withdraw()
        BackButton = Button(displayChain, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back_1, bg="red", fg="white", relief=RAISED)
    elif Option == 2:
        gene = Protein
        addOwn.withdraw()
        BackButton = Button(displayChain, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back_2, bg="red", fg="white", relief=RAISED)
        
    chain = gene_translation(gene)[1] #O(n**2)
    BackButton.place(relx=0.1, rely=0.1, anchor=CENTER)

    # Define grid dimensions
    rows = 9
    cols = 9
    cell_width = 55
    cell_height = 55 
    padding_x = 10
    padding_y = 10

    # Draw circles and place text inside each circle
    if len(chain) == 0:
        not_found = Label(displayChain, text="No amino acids found...", font=("Times New Roman", 20, "bold"), bg= "white", fg="red", relief=FLAT)
        not_found.place(anchor=CENTER, relx=0.5, rely=0.5)
    else:
        try:
            for i in range(rows): #O(n**2)
                for j in range(cols):
                    index = i * cols + j
                    x = padding_x + j * (cell_width + padding_x * 2)
                    y = padding_y + i * (cell_height + padding_y * 2)
                    draw_circle_with_text(canvas, x + 30, y + 30, 30, chain[index]) #O(1)
        except:
            pass
    
    proceed_button = Button(displayChain, text="Display Sequence", font=("Times New Roman", 13, "bold"), width=23, height=2, command=lambda: open_Display_List(Protein,Option), bg="purple", fg="white", relief=SUNKEN)
    proceed_button.place(relx=0.98, rely=0.85, anchor=SE)

    add_pattern = Button(displayChain, text="Enter Pattern To Search", font=("Times New Roman", 13, "bold"), width=23, height=2, command=lambda: Input_Pattern(Protein,Option), bg="white", fg="purple", relief=SUNKEN)
    add_pattern.place(anchor=SE, relx = 0.98, rely = 0.95)
 
def open_Display_List(Protein, Option): #display each base with its relative index in string; O(n**2)

    global displayStrand, displayChain, bg_image, compare_root

    displayStrand = Toplevel()
    displayStrand.title("Gene Sequence Display")
    displayStrand.geometry("1920x1080")

    if Option == 1:
        gene_sequence = csv_reader(Protein) #O(n+m)
    else:
        gene_sequence = Protein

    bg_label = Label(displayStrand, image=bg_image)
    bg_label.place(relx=0, rely=0, relwidth=1, relheight=1)

    # Create a label to display the gene name if we know it
    if Option == 1:
        gene_name_label = Label(displayStrand, text=Protein, font=("Times New Roman", 14, "bold"), width = 10 + 2*len(Protein), fg="white", bg="black", pady=10, relief=SUNKEN)
        gene_name_label.place(anchor=CENTER, relx=0.5, rely=0.1)
    else:
        gene_name_label = Label(displayStrand, text="Your Sequence: ", font=("Times New Roman", 14, "bold"), width = 40, fg="white", bg="black", pady=10, relief=SUNKEN)
        gene_name_label.place(anchor=CENTER, relx=0.5, rely=0.1)

    length = len(gene_sequence) 

    # Add empty bases so that we may print full lines of 21 bases
    if length < 252: 
        for _ in range(252 - length): # 243 chars max so we use 252 to make step of 21 for the loop
            gene_sequence += "-" 

    for i in range(0, 252, 21): #O(n**2)
        line = gene_sequence[i:i+21]
        if line == "-"*21: #break when entire line is empty
            break
        index_line = ''.join(f"{line[j]}{[i+j]} "for j in range(21) if line[j] != "-")
        sequence_label = Label(displayStrand, text=index_line, font=("Courier", 10, "bold"), fg="black")
        sequence_label.place(relx=0.5, rely=0.2+(0.05*(i//21)), anchor=CENTER)

    try:
        displayChain.withdraw()
        compare_root.withdraw()
    except:
        pass

    def back(): #return to polypeptide chain
        displayStrand.destroy()
        displayChain.deiconify()
    
    def end(): #finish searching for pattern
        displayStrand.destroy()
        searchroot.destroy()

    def end2(): #finish after mutating (Final)
        displayStrand.destroy()
        compare_root.destroy()

    if Option == 1 or Option == 2:
        TopButton = Button(displayStrand, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back, bg="red", fg="white", relief=RAISED)
    elif Option == 3:
        TopButton = Button(displayStrand, text="END", font=("Times New Roman", 16), width=10, height=1, command=end, bg="red", fg="white", relief=RAISED)
    elif Option == 4:
        TopButton = Button(displayStrand, text="END", font=("Times New Roman", 16), width=10, height=1, command=end2, bg="red", fg="white", relief=RAISED)

    TopButton.place(relx=0.1, rely=0.1, anchor=CENTER)

def Input_Pattern(Protein, Option): #input pattern to search for; O(n**2)

    global editStrand, displayChain, bg_image, compare_root, sequence, bwt, suffix_array, c_table, occ_ranks, mutation_count, first_sequence, coordinates

    editStrand = Toplevel()
    editStrand.title("Input Pattern")
    editStrand.geometry("1920x1080")

    if Option == 1:
        gene_sequence = csv_reader(Protein).upper() #O(n+m)
    else:
        gene_sequence = Protein.upper()

    if mutation_count == 0:
        first_sequence = gene_sequence.upper()

    sequence = gene_sequence.upper()

    bg_label = Label(editStrand, image=bg_image)
    bg_label.place(relx=0, rely=0, relwidth=1, relheight=1)

    input_frame = Frame(editStrand)
    input_frame.place(anchor=CENTER, relx=0.5, rely=0.5)
    
    def process_DisplayList(pattern):
        global sequence, bwt, suffix_array, c_table, occ_ranks, coordinates
        pattern = pattern.strip()
        if len(pattern) < 3:
            instruct_Label.config(text="Pattern should be atleast 3 bases long!", fg="red")
        elif len(pattern) > len(sequence):
            instruct_Label.config(text="Pattern longer than Sequence", fg="red")
        else:
            result = validation_check(pattern) #O(n)
            if result == "Processing...":
                bwt, suffix_array = bwt_and_suffix_array(sequence) #O(n)
                c_table = create_c_table(bwt) #O(nlogn)
                occ_ranks = rank_bwt(bwt) #O(n**2)
                coordinates = search_with_fm_index(pattern, bwt, c_table, occ_ranks, suffix_array) #O(nlogn)
                if coordinates == -1:
                    instruct_Label.config(text="Pattern not found...", fg="red")
                else:
                    instruct_Label.config(text="Enter a template pattern to search:", fg="white")
                    Search()
            else:
                instruct_Label.config(text=result, fg="red")
    
    try:
        displayChain.withdraw()
        compare_root.withdraw()
    except:
        pass

    # Label to display the output message
    instruct_Label = Label(editStrand, text="Enter a template pattern to search:", font=("Times New Roman", 14, "bold"), width=40, height=1, bg='black',fg='white')
    instruct_Label.place(anchor= CENTER, relx=0.5, rely = 0.3)
    
    # Entry widgets for entering starting indexes
    pattern = Entry(input_frame, width=30)
    pattern.grid(row=0, column=0, padx=2)

    # Button to process the input
    process_button = Button(input_frame, text="Search", command=lambda:process_DisplayList(pattern.get().upper()), bg="black", fg="white")
    process_button.grid(row=0, column=9, padx=6)

    def back1():
        editStrand.withdraw()
        displayChain.deiconify()

    def back2():
        editStrand.withdraw()
        compare_root.deiconify()

    if Option == 3:
        BackButton = Button(editStrand, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back2, bg="red", fg="white", relief=RAISED)
    else:
        BackButton = Button(editStrand, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back1, bg="red", fg="white", relief=RAISED)
    BackButton.place(relx=0.1, rely=0.1, anchor=CENTER)

def Search(): #Returns indices where pattern was found; O(n**2)

    global searchroot, bg_image, editStrand, sequence, coordinates

    searchroot = Toplevel()
    searchroot.geometry("1920x1080")
    searchroot.title("Pattern Search")

    bg_label = Label(searchroot, image=bg_image)
    bg_label.place(relx=0, rely=0, relwidth=1, relheight=1)
    
    heading = Label(searchroot, text = "Your pattern was found at the following indices: ", font=("Times New Roman", 16, "bold"), width=50, height=1, bg="pink", fg="black", relief=FLAT)
    heading.place(relx = 0.5, rely = 0.1, anchor=CENTER)

    length = len(coordinates)

    # Add placeholders "-" to make sure the total length is a multiple of 11
    for _ in range(253 - length):
        coordinates += "-"

    # Calculate the number of lines needed
    num_lines = (length + 12) // 11

    # Display tuples in lines of 11 tuples each
    for i in range(num_lines): #overall O(11*11*n) => O(n)
        line = coordinates[i * 11: (i + 1) * 11]  # Extract tuples for this row
        if all(tup == "-" for tup in line): #O(11)
            break
        index_line = ' '.join(str(tup) for tup in line if tup != "-") #O(11)
        sequence_label = Label(searchroot, text=index_line, font=("Courier", 12, "bold"), fg="black")
        sequence_label.place(relx=0.5, rely=0.15 + (0.03 * i), anchor=CENTER)

    End = Button(searchroot, text = "Proceed To End", command=lambda: open_Display_List(sequence, 3), font=("Times New Roman", 14), width=20, height=2, bg="bisque", fg="black", relief=RAISED)
    Mut = Button(searchroot, text = "Proceed To Mutations", command=lambda: Mutate(),  font=("Times New Roman", 14), width=20, height=2, bg="yellow", fg="black", relief=RAISED)
    End.place(relx = 0.3, rely = 0.85, anchor=CENTER)
    Mut.place(relx = 0.7, rely = 0.85, anchor=CENTER)
    
    editStrand.withdraw()

    def back():
        searchroot.withdraw()
        editStrand.deiconify()

    BackButton = Button(searchroot, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back, bg="red", fg="white", relief=RAISED)
    BackButton.place(relx=0.1, rely=0.1, anchor=CENTER)

def Mutate(): #Takes indices as input and performs mutations there; O(n**2)

    global mutation_root, bg_image, searchroot, coordinates

    mutation_root = Toplevel()
    mutation_root.geometry("1920x1080")
    mutation_root.title("Time To Mutate!")

    bg_label = Label(mutation_root, image=bg_image)
    bg_label.place(relx=0, rely=0, relwidth=1, relheight=1)

    searchroot.withdraw()

    def validate_coords(coords, coordinates): #Ensure tuple is entered and it exists in list

        # Assume entry is a tuple
        tup_bool = True
        # Assume entered tuple exists
        found_bool = True

        #destroy previous error messages when new entry is made
        for widget in mutation_root.winfo_children(): #O(n)
            if isinstance(widget, Label) and widget.cget("bg") == "black" and widget.cget("fg") == "red":
                widget.destroy()
            
        try:
            a, b = map(int, coords.strip().split(",")) #O(n)
            coordtuple = tuple((a, b))  # Assign coords as a tuple of integers
        except ValueError:
            display_error("Please enter two integers separated by a comma!")
            tup_bool = False

        # If the entry is a valid tuple, check if it exists in the coordinates
        if tup_bool:
            if coordtuple not in coordinates:
                found_bool = False
                # Display error label if tuple does not exist in coordinates
                display_error("Tuple Does Not Exist!")

        # If both tuple format and existence are correct, proceed with mutation process
        if tup_bool and found_bool:
            val_coords.destroy() #remove validate button
            coords_entry.destroy() #remove entry box

            add_button = Button(inputmut_frame, text="Insertion", command=lambda: process_mutate(coordtuple, "Insertion"), bg="yellow", fg="black", width=12, height=2)
            add_button.grid(row=1, column=0, padx=10, pady=10)

            sub_button = Button(inputmut_frame, text="Substitution", command=lambda: process_mutate(coordtuple, "Substitution"), bg="green", fg="black", width=12, height=2)
            sub_button.grid(row=1, column=1, padx=10, pady=10)

            del_button = Button(inputmut_frame, text="Deletion", command=lambda: process_mutate(coordtuple, "Deletion"), bg="cyan", fg="black", width=12, height=2)
            del_button.grid(row=1, column=2, padx=10, pady=10)

    def display_error(message):
        error = Label(mutation_root, text=message, bg="black", fg="red", width=10 + len(message), height=1)
        error.place(anchor=CENTER, relx=0.5, rely=y+0.15)

    def process_mutate(coordtuple, command):

        global sequence, bwt, suffix_array, c_table, occ_ranks

        for widget in inputmut_frame.winfo_children(): #O(n)
            if isinstance(widget, Button):
                if widget.cget("text") in ["Insertion", "Substitution", "Deletion"]:
                    widget.destroy()

        og_sequence = sequence #save pre-mutation sequence

        if command == 'Insertion' or command == 'Substitution':
            subSequence = Entry(inputmut_frame, width=30)
            subSequence.grid(row=1, column=1, padx=10, pady=10)
            val_subSeq = Button(inputmut_frame, text = "Proceed "+command, width=20, command=lambda:process_sub(subSequence.get(), coordtuple, command, og_sequence), bg = "black", fg="white")
            val_subSeq.grid(row=2, column=1, padx=10, pady=10)
        else:
            #if deletion, no need to take subsequence entry
            sequence, bwt, suffix_array, c_table, occ_ranks = mutations(sequence, "Deletion", coordtuple) #O(n**2)
            compare(og_sequence, sequence)

    def process_sub(subSequence, coordtuple, command, og_sequence):

        global sequence, bwt, suffix_array, c_table, occ_ranks

        if command == "Insertion" and len(subSequence) + len(sequence) > 243:
            error = Label(mutation_root, text = "Sequence can not exceed 243 bases!", bg= "black", fg = "red", width= 30, height=1)
            error.place(anchor=CENTER, relx=0.5, rely = y+0.17)
        else:
            for widget in mutation_root.winfo_children(): #O(n)
                if isinstance(widget, Label) and widget.cget("bg") == "black" and widget.cget("fg") == "red":
                    widget.destroy()

            if validation_check(subSequence) == "Processing...":
                sequence, bwt, suffix_array, c_table, occ_ranks = mutations(sequence, command, coordtuple, subSequence.upper()) #O(n**2)
                compare(og_sequence, sequence)
            else:
                error = Label(mutation_root, text = "Invalid SubSequence Entry!", bg= "black", fg = "red", width= 25, height=2)
                error.place(anchor=CENTER, relx=0.5, rely = y+0.17)

    given_coords = Label(mutation_root, text='Please enter one of the tuples pattern was found at: ', font=("Times New Roman", 16), width=50, height=1, bg="pink", fg="black", relief=FLAT)
    given_coords.place(relx=0.5, rely=0.1, anchor=CENTER)

    length = len(coordinates)

    # Add placeholders "-" to make sure the total length is a multiple of 11
    for _ in range(253 - length):
        coordinates += "-"

    # Calculate the number of lines needed
    num_lines = (length + 12) // 11

    # Display tuples in lines of 11 tuples each
    for i in range(num_lines): #overall O(11*11*n) => O(n)
        line = coordinates[i * 11: (i + 1) * 11]  # Extract tuples for this row
        if all(tup == "-" for tup in line): #O(11)
            break
        index_line = ' '.join(str(tup) for tup in line if tup != "-") #O(11)
        sequence_label = Label(mutation_root, text=index_line, font=("Courier", 12, "bold"), fg="black")
        y = 0.15 + (0.03 * i)
        sequence_label.place(relx=0.5, rely=0.15 + (0.03 * i), anchor=CENTER)

    inputmut_frame = Frame(mutation_root)
    inputmut_frame.place(anchor=CENTER, relx=0.5, rely=y+0.08)

    coords_entry = Entry(inputmut_frame, width=30)
    coords_entry.grid(row=0, column=0, padx=10, pady=10)

    val_coords = Button(inputmut_frame, text = "Mutate", width=15, command=lambda:validate_coords(coords_entry.get(), coordinates), bg = "black", fg="white")
    val_coords.grid(row = 0, column = 2, padx=10, pady=10)

    def back():
        mutation_root.destroy()
        searchroot.deiconify()

    back_button = Button(mutation_root, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back, bg="red", fg="white", relief=RAISED)
    back_button.place(relx=0.1, rely=0.1, anchor=CENTER)

def compare (og, new): #Draw two canvases to show chain before and after mutation; O(n**2)

    global compare_root, mutation_root, bg_image, searchroot, mutation_count, first_sequence

    compare_root = Toplevel()
    compare_root.title("Compare PolyPeptides!")
    compare_root.geometry("1920x1080")

    bg_label = Label(compare_root, image=bg_image)
    bg_label.place(relx=0, rely=0, relwidth=1, relheight=1)

    mutation_root.withdraw()

    def back():
        compare_root.withdraw()
        searchroot.deiconify()

    mutation_count += 1
    
    if mutation_count == 3:
        og = first_sequence
        prevlabel = "Original Sequence:"
    else:
        prevlabel = "Previous Sequence:"

    back_button = Button(compare_root, text="Back", font=("Times New Roman", 16), width=10, height=1, command=back, bg="red", fg="white", relief=RAISED)
    back_button.place(relx=0.1, rely=0.1, anchor=CENTER)

    canvas_og = Canvas(compare_root, bg="white") 
    canvas_og.place(relx=0.25, rely=0.5, relwidth=0.45, relheight=0.6, anchor=CENTER)

    canvas_new = Canvas(compare_root, bg="white") 
    canvas_new.place(relx=0.75, rely=0.5, relwidth=0.45, relheight=0.6, anchor=CENTER)

    mutation_no = Label(compare_root,text = "Mutation number " + str(mutation_count), font=("Times New Roman", 16, "bold", "underline"), width=20, height=1, bg="pink", fg="black", relief=FLAT)
    mutation_no.place(relx=0.5, rely=0.08, anchor = CENTER)

    #Prepare error message that stop codon has been produced so the full chain will not be shown
    polypeptidenew = gene_translation(new) #O(n**2)
    if len(polypeptidenew) == 3:
        nonsense_msg = polypeptidenew[2]
    else:
        nonsense_msg = None
    og_chain = gene_translation(og)[1] #O(n**2)
    new_chain = polypeptidenew[1]

    rows = 9 
    cols = 11
    cell_width = 50
    cell_height = 50
    padding_x = 2
    padding_y = 2

    og_label = Label(compare_root, text = prevlabel, font=("Times New Roman", 16, "bold", "underline"), width=20, height=1, bg="purple", fg="white", relief=FLAT)
    new_label = Label(compare_root, text = "Mutated Sequence:", font=("Times New Roman", 16, "bold", "underline"), width=20, height=1, fg="purple", bg="white", relief=FLAT)
    og_label.place(anchor=CENTER, relx= 0.25, rely= 0.15)
    new_label.place(anchor=CENTER, relx= 0.75, rely= 0.15)

    if nonsense_msg is not None:
        mutation_label = Label(compare_root, text=nonsense_msg, font=("Times New Roman", 16, "bold", "underline"), width=50, height=1, bg="white", fg="black", relief=FLAT)
        mutation_label.place(anchor=CENTER, relx= 0.75, rely= 0.85)

    #display original chain; O(n**2)
    try:
        for i in range(rows):
            for j in range(cols):
                index = i * cols + j
                x_og = padding_x + j * (cell_width + padding_x * 2)
                y_og = padding_y + i * (cell_height + padding_y * 2)
                draw_circle_with_text(canvas_og, x_og + 10, y_og + 10, 10, og_chain[index])
    except:
        pass

    #display new chain; O(n**2)
    #new loop is needed for this because lengths of og and new sequence are not same
    try:
        for i in range(rows):
            for j in range(cols):
                index = i * cols + j
                x_new = padding_x + j * (cell_width + padding_x * 2)
                y_new = padding_y + i * (cell_height + padding_y * 2)
                draw_circle_with_text(canvas_new, x_new + 10, y_new + 10, 10, new_chain[index])
    except:
        pass

    # make sure user can only make at most 3 mutations to a single protein
    def process_again():

        global mutation_count

        if mutation_count < 3:
            Input_Pattern(new, 3)
        else:
            open_Display_List(new, 4)

    if mutation_count == 3:
        again_text = "End Sequence"
    else:
        again_text = "Make Another Mutation"

    again = Button(compare_root, text= again_text, font=("Times New Roman", 14), width=20, height=1, command=lambda:process_again(), bg="white", fg="black", relief=RAISED)
    end = Button(compare_root, text= "End Program", font=("Times New Roman", 14), width=20, height=1, command=lambda: open_Display_List(new, 4), fg="white", bg="black", relief=RAISED)
    again.place(anchor=CENTER, relx=0.4, rely= 0.9)
    end.place(anchor=CENTER, relx=0.6, rely = 0.9)

main()

# Overall Time Complexity O(n**2)
# Overall Space Complexity O(n)
# Learning Reference: "Python GUI: Tkinter Tutorial In Hindi" by CodeWithHarry (https://youtube.com/playlist?list=PLu0W_9lII9ajLcqRcj4PoEihkukF_OTzA&si=AneQcjB1Gm2QAwDk)