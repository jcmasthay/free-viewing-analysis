function t = convert_char_vars_to_string(t)

for i = 1:size(t, 2)
  tvar = t.(t.Properties.VariableNames{i});
  t.(t.Properties.VariableNames{i}) = convertCharsToStrings( tvar );
end

end