# TODO write documentation

function specialchar_cancel_line()
    return "\r\u001b[2K"
end

return specialchar_cancel_line

function specialchar_move_cursor_left_by_1()
    return "\u001b[1D"
end

return specialchar_move_cursor_left_by_1

function specialchar_backspace()
    return "\b"
end

return specialchar_backspace

function specialchar_newline()
    print("\n")
end

return specialchar_newline