function s = session_to_datetime(s)

s = datetime( datestr(datenum(s, 'mmddyyyy')) );

end