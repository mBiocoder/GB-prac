package com.company;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class Parser {

    private HashMap<String, String> parameterMap;
    private HashSet<String> allowedOptions;
    private HashSet<String> obligatoryOptions;
    private HashSet<String> nonParameterOptions;

    private HashSet<String> intOptions;
    private HashSet<String> doubleOptions;
    private HashSet<String> fileOptions;

    //Constructor
    public Parser(String... allowedOptions) {
        this.allowedOptions = new HashSet<>();
        this.allowedOptions.addAll(Arrays.asList(allowedOptions));
        this.obligatoryOptions = new HashSet<>();
        this.parameterMap = new HashMap<>();
        this.intOptions = new HashSet<>();
        this.doubleOptions = new HashSet<>();
        this.fileOptions = new HashSet<>();
        this.nonParameterOptions = new HashSet<>();
    }

    public void setObligatoryOptions(String... options) {
        this.obligatoryOptions.addAll(Arrays.asList(options));
    }

    public void parseArguments(String[] args) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith("-")) {
                if (!allowedOptions.contains(args[i])) {
                    throw new RuntimeException("Parameter not allowed: " + args[i]);
                } else if (i + 1 == args.length && !nonParameterOptions.contains(args[i])) {
                    throw new RuntimeException("Last parameter has no value: " + args[i]);
                } else if (parameterMap.containsKey(args[i])) {
                    throw new RuntimeException("Parameter already in map: " + args[i]);
                    // check if argument with '-' is a new argument or a negative number
                } else if (i == args.length-1 || args[i + 1].startsWith("-") && !Arrays.asList(new String[]{"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"}).contains(args[i + 1].substring(1, 2))) {
                    if (nonParameterOptions.contains(args[i])) {
                        parameterMap.put(args[i], null);
                    } else {
                        throw new RuntimeException("Parameter without value: " + args[i]);
                    }
                } else {
                    parameterMap.put(args[i], args[++i]);
                }
            } else {
                throw new RuntimeException("Argumets have to start with '-'");
            }
        }

        // check if all obligatory options are there
        for (String option : obligatoryOptions) {
            if (!parameterMap.containsKey(option)) {
                throw new RuntimeException("Obligatory option missing: " + option);
            }
        }
    }

    // set parameter types
    public void setInt(String... intOptions) {
        this.intOptions.addAll(Arrays.asList(intOptions));
    }

    public void setDouble(String... doubleOptions) {
        this.doubleOptions.addAll(Arrays.asList(doubleOptions));
    }

    public void setFile(String... fileOptions) {
        this.fileOptions.addAll(Arrays.asList(fileOptions));
    }

    public void setNonParameterOptions(String... nonParameterOptions) {
        this.nonParameterOptions.addAll(Arrays.asList(nonParameterOptions));
    }

    public String getValue(String option) {
        return parameterMap.getOrDefault(option, null);
    }

    public Integer getInt(String option) {
        if (!parameterMap.containsKey(option)) {
            return null;
        } else if (intOptions.contains(option)) {
            return Integer.parseInt(getValue(option));
        } else {
            throw new RuntimeException("Option value is not an Integer: " + option);
        }
    }

    public Double getDouble(String option) {
        if (!parameterMap.containsKey(option)) {
            return null;
        } else if (doubleOptions.contains(option)) {
            return Double.parseDouble(getValue(option));
        } else {
            throw new RuntimeException("Option value is not a Double: " + option);
        }
    }

    public File getFile(String option) {
        if (!parameterMap.containsKey(option)) {
            return null;
        } else if (fileOptions.contains(option)) {
            return new File(getValue(option));
        } else {
            throw new RuntimeException("Option value is not a File: " + option);
        }
    }

    public boolean isSet(String option) {
        if (!parameterMap.containsKey(option)) {
            return false;
        } else if (nonParameterOptions.contains(option)) {
            return true;
        } else {
            throw new RuntimeException("Option value is not a File: " + option);
        }
    }

}