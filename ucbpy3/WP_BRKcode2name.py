def code2name( code_list = ['-11','-12']):
    """
    SW naming convention taken from
    InversionCodes/Libs/ucblibs/include/structwavket.h
    """
    labels_name = []
    for lb in code_list:
        phase_name = ''
        # check that the wp corresponds to SW
        if lb[0] == '-':
            assert len(lb) == 3, "WP identifier should have 3 characters"
            if float(lb[1]) < 5:
                if lb[1] == '1':
                     phase_name += 'R'
                elif lb[1] == '2':
                     phase_name += 'G'
                elif lb[1] == '3':
                     phase_name += 'XR'
                elif lb[1] == '4':
                     phase_name += 'XG'

                phase_name += lb[2]
            else:
                if lb[2] == '1':
                    p1 = '1'
                    p2 = '1'
                elif lb[2] == '2':
                    p1 = '1'
                    p2 = '2'
                elif lb[2] == '3':
                    p1 = '2'
                    p2 = '2'
                elif lb[2] == '4':
                    p1 = '2'
                    p2 = '3'
                elif lb[2] == '5':
                    p1 = '3'
                    p2 = '4'
 
                if lb[1] == '5':
                    phase_name = 'R'+p1+'+X'+p2
                elif lb[1] == '6':
                    phase_name = 'G'+p1+'+GX'+p2
  
                if lb[1:] == '71':
                    phase_name = 'RX1+RX2'
                if lb[1:] == '72':
                    phase_name = 'RX2+RX3'

                if lb[1:] == '81':
                    phase_name = 'GX1+GX2'
                if lb[1:] == '82':
                    phase_name = 'GX2+GX3'

        else:
           if lb[:] == "20":
               phase_name = "Sdiff"
           elif lb[:] == "21":
               phase_name = "S"
           elif lb[:] == "22":
               phase_name = "SS"
           elif lb[:] == "23":
               phase_name = "SSS"
           elif lb[:] == "24":
               phase_name = "S4"
           elif lb[:] == "25":
               phase_name = "ScS"
           elif lb[:] == "26":
               phase_name = "ScS$_2$"
           elif lb[:] == "27":
               phase_name = "ScS$_3$"
           elif lb[:] == "44":
               phase_name = "sS$_3$"
           elif lb[:] == "121":
               phase_name = "S+ScS"
           elif lb[:] == "123":
               phase_name = "SSS+ScS$_2$"
           elif lb[:] == "124":
               phase_name = "sS+sScS"
           elif lb[:] == "126":
               phase_name = "SSS+sSS"
           elif lb[:] == "127":
               phase_name = "S$_4$+sS$_3$"
           elif lb[:] == "144":
               phase_name = "ScS+SS"
           elif lb[:] == "148":
               phase_name = "S$_4$+ScS$_2$"
           elif lb[:] == "1000":
               phase_name = "unknown"
           # From the BW picking code
           elif lb[:] == "100":
               phase_name = "S"
           elif lb[:] == "100":
               phase_name = "Sdiff"
           elif lb[:] == "101":
               phase_name = "S"
           elif lb[:] == "102":
               phase_name = "SS"
           elif lb[:] == "103":
               phase_name = "P"
           elif lb[:] == "104":
               phase_name = "Pdiff"
           elif lb[:] == "105":
               phase_name = "SSS"
           elif lb[:] == "106":
               phase_name = "S$_4$"
           elif lb[:] == "800":
               phase_name = "unknown"



           else:
               phase_name = lb[:]






        labels_name.append( phase_name)
    return labels_name


