#!/usr/bin/env python
# -*- coding: utf-8 -*-

# see. http://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output

def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")

    return message

                    
class ColoredFormatter(logging.Formatter):
    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[1;{}m"
    BOLD_SEQ = "\033[1m"

    COLORS = {
        'WARNING': YELLOW,
        'INFO': WHITE,
        'DEBUG': BLUE,
        'CRITICAL': YELLOW,
        'ERROR': RED
    }

    def __init__(self, msg, use_color =True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if (self.use_color) and (levelname in ColoredFormatter.COLORS):
            levelname_color = ColoredFormatter.COLOR_SEQ.format(30 + ColoredFormatter.COLORS[levelname]) + levelname + ColoredFormatter.RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)

