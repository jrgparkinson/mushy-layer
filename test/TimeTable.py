import os
import re

class TimeTableMethod:

    def __init__(self, string_to_parse, parent=None):
        parts = re.findall('(\s*)\[(\d+)]\s([^\s]*)\s([\d\.]+)\s+[\d\.]+\%\s', string_to_parse)

        self.parent = parent
        self.valid = False

        if parts:
            match = parts[0]
            # print(match)
            self.indent = len(match[0])
            self.id = int(match[1])
            self.name = match[2].strip()
            self.time = float(match[3])
            self.valid = True
        # else:
        #     print('Could not parse: %s' % string_to_parse)


    def __repr__(self):
        # return '%s/%s' % (self.parent, self.name)
        return self.name

    def long_desc(self):
        return '[%d] %s (%s)' % (self.id, self.name, self.time)

class TimeTable:

    def __init__(self, filepath):
        self.filepath = filepath

        self.methods = []

        self.parse_file()



    def parse_file(self):

        if os.path.exists(self.filepath):
            with open(self.filepath, 'r') as f:

                file_contents = f.read()

                # print(file_contents)

                # Split into methods
                # matches = re.findall('[\-]+\n([^\-])', file_contents)

                # Find the second [0] tag
                parts = file_contents.split('[0]')
                second_half = parts[2]
                lines = second_half.split('\n')
                # print(lines)


                main_method = TimeTableMethod('[0]' + lines[0])
                current_parent = main_method

                prev_method = main_method

                for line in lines[1:]:
                    method = TimeTableMethod(line)

                    if not method.valid:
                        continue

                    # If this method is more indented, change parent
                    if method.indent > prev_method.indent:
                        current_parent = prev_method

                    # If this method is less indented, parent increases
                    if method.indent < prev_method.indent:
                        indent_change = (prev_method.indent - method.indent)
                        while (indent_change > 0) and current_parent.parent:
                            current_parent = current_parent.parent

                            # Indent changes by 3 per level
                            indent_change = indent_change - 3

                    method.parent = current_parent

                    self.methods.append(method)

                    prev_method = method


                # for method in self.methods[0:5]:
                #     print('%s' % method)


                # for line in f.readlines():
                #
                #     # New methods start with '['
                #     if line[0] == '-':


    def get_all_children(self, parent_method_name):
        '''
        :param parent_method_name:
        :return: dict of individual parent occurences mapping to their children
        '''


        parent_children = {}

        for m in self.methods:
            if m.parent and m.parent.name == parent_method_name:
                # children.append(m)
                if m.parent not in list(parent_children.keys()):
                    parent_children[m.parent] = []

                parent_children[m.parent].append(m)


        return parent_children

    def total_time_in_method(self, method_name):
        total_time = 0
        for m in self.methods:
            if m.name == method_name:
                total_time = total_time + m.time

        return total_time

    def get_call_history_for_id(self, method_id):
        call_history = []

        for m in self.methods:
            if m.id == method_id:
                call_history.append(m)

                p = m.parent
                call_history.append(p)
                while p.parent:
                    p = p.parent
                    call_history.append(p)

                call_history = call_history[::-1]
                break

        return call_history


    def get_call_history_for_name(self, method_name):

        call_history = []

        for m in self.methods:
            if m.name == method_name:
                call_history.append(self.get_call_history_for_id(m.id))

        return call_history




if __name__ == "__main__":
    filepath = 'time.table.0.example'
    print('Example usage with file: %s' % filepath)

    time_table = TimeTable(filepath)

    method_name = 'VCAMRPoissonOp2::restrictResidual'
    method_name = 'AMRLevelMushyLayer::calculateTimeIndAdvectionVel'
    # method_name = 'AMR::timeStep'
    method_name = 'AMRNonLinearMultiCompOp::levelGSRB::BCs'
    parent_children = time_table.get_all_children(method_name)

    print('Children of %s\n' % method_name)
    for parent in list(parent_children.keys()):
        print(parent.long_desc())
        for c in parent_children[parent]:
            print('  %s' % c.long_desc())

    print('Total time in %s: %.3g' % (method_name, time_table.total_time_in_method(method_name)))

    print('Call histories for %s:' % method_name)
    call_histories = time_table.get_call_history_for_name(method_name)

    for c in call_histories:

        str = ' -> '.join(['%s' % m for m in c])
        print(str + '\n')



    # Compare multiple files

    files = []
    tt = [TimeTable(f) for f in files]

    


