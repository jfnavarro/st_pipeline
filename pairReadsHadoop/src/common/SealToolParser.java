// Copyright (C) 2011-2012 CRS4.
//
// This file is part of Seal.
//
// Seal is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Seal is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License along
// with Seal.  If not, see <http://www.gnu.org/licenses/>.

package common;

import common.ClusterUtils;

import java.io.IOException;
import java.io.File;
import java.io.FileReader;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.cli.*;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.util.GenericOptionsParser;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public class SealToolParser {

	public static final File DefaultConfigFile = new File(System.getProperty("user.home"), ".sealrc");
	public static final int DEFAULT_MIN_REDUCE_TASKS = 0;
	public static final int DEFAULT_REDUCE_TASKS_PER_NODE = 3;

	public static final String INPUT_FORMAT_CONF = "seal.input-format";
	public static final String OUTPUT_FORMAT_CONF = "seal.output-format";

	protected static final String INPUT_FORMAT_DESC = "Input format name";
	protected static final String OUTPUT_FORMAT_DESC = "Output format name";

	private int minReduceTasks;

	/**
	 * Configuration object used to parse the command line, cached for further queries.
	 */
	private Configuration myconf;

	protected Options options;
	private Option opt_nReduceTasks;
	private Option opt_configFileOverride;
	private Option opt_inputFormat;
	private Option opt_outputFormat;
	private Option opt_compressOutput;

	private String[] acceptedInputFormats;
	private String[] acceptedOutputFormats;

	private Integer nReduceTasks;
	private int nReduceTasksPerNode;
	private String configSection;
	protected String toolName;

	protected ArrayList<Path> inputs;
	private Path outputDir;

	/**
	 * Construct a SealToolParser instance.
	 *
	 * The instance is set to read the properties in configuration file's section sectionName,
	 * in addition to the default section.  Properties set on the command line will override
	 * the file's settings.
	 *
	 * @param configSection Name of section of configuration to load, in addition to DEFAULT.
	 * If null, only DEFAULT is loaded
	 * @param toolName Name used in the help message
	 */
	public SealToolParser(String configSection, String toolName)
	{
		this.toolName = toolName;

		options = new Options(); // empty
		opt_nReduceTasks = OptionBuilder
			              .withDescription("Number of reduce tasks to use.")
			              .hasArg()
			              .withArgName("INT")
			              .withLongOpt("num-reducers")
			              .create("r");
		options.addOption(opt_nReduceTasks);

		opt_configFileOverride = OptionBuilder
			              .withDescription("Override default Seal config file (" + DefaultConfigFile + ")")
			              .hasArg()
			              .withArgName("FILE")
			              .withLongOpt("seal-config")
			              .create("sc");
		options.addOption(opt_configFileOverride);

		opt_inputFormat = OptionBuilder
			              .withDescription(INPUT_FORMAT_DESC)
			              .hasArg()
			              .withArgName("FORMAT")
			              .withLongOpt("input-format")
			              .create("if");
		options.addOption(opt_inputFormat);

		opt_outputFormat = OptionBuilder
			              .withDescription(OUTPUT_FORMAT_DESC)
			              .hasArg()
			              .withArgName("FORMAT")
			              .withLongOpt("output-format")
			              .create("of");
		options.addOption(opt_outputFormat);

		opt_compressOutput = OptionBuilder
			              .withDescription("Compress output files with CODEC (one of gzip, bzip2, snappy, auto)")
			              .hasArg()
			              .withArgName("CODEC")
			              .withLongOpt("compress-output")
			              .create("oc");
		options.addOption(opt_compressOutput);

		nReduceTasks = null;
		inputs = new ArrayList<Path>(10);
		outputDir = null;
		this.configSection = (configSection == null) ? "" : configSection;
		minReduceTasks = DEFAULT_MIN_REDUCE_TASKS;
		myconf = null;
		nReduceTasksPerNode = DEFAULT_REDUCE_TASKS_PER_NODE;
	}

	/**
	 * Set the minimum acceptable number of reduce tasks.
	 * If a user specifies a number lower than this limit parseOptions will raise
	 * an error.
	 */
	public void setMinReduceTasks(int x)
	{
		if (x < 0)
			throw new IllegalArgumentException("minimum number of reduce tasks must be >= 0");
		minReduceTasks = x;
	}

	public int getMinReduceTasks() { return minReduceTasks; }

	protected void loadConfig(Configuration conf, File fname) throws ParseException, IOException
	{
		ConfigFileParser parser = new ConfigFileParser();

		try {
			parser.load( new FileReader(fname) );

			Iterator<ConfigFileParser.KvPair> it = parser.getSectionIterator(configSection);
			ConfigFileParser.KvPair pair;
			while (it.hasNext())
			{
				pair = it.next();
				conf.set(pair.getKey(), pair.getValue());
			}
		}
		catch (FormatException e)
		{
			throw new ParseException("Error reading config file " + fname + ". " + e);
		}
	}

	/**
	 * Decides whether to use an rc file, and if so which one.
	 *
	 * This method is necessary only because we'd like the user to be able to override the default
	 * location of the seal configuration file ($HOME/.sealrc).  So, it scans
	 * the command line arguments looking for a user-specified seal configuration file.
	 * If one is specified, it verifies that it exists and is readable.  If none is specified
	 * it checks to see whether a configuration file is available at the default location,
	 * and if it is the method verifies that it is readable.
	 *
	 * If a config file is found and is readable, its path is returned as a File object.  On the other
	 * hand, if a config file isn't found the method returns null.
	 *
	 * @param args command line arguments
	 * @exception ParseException raise if the file specified on the cmd line doesn't exist or isn't readable.
	 */
	protected File getRcFile(String[] args) throws ParseException
	{
		File fname = null;

		String shortOpt = "--" + opt_configFileOverride.getOpt();
		String longOpt = "--" + opt_configFileOverride.getLongOpt();

		for (int i = 0; i < args.length; ++i)
		{
			if (args[i].equals(shortOpt) || args[i].equals(longOpt))
			{
				if (i+1 >= args.length)
					throw new ParseException("Missing file argument to " + args[i]);

				fname = new File(args[i+1]);
				break;
			}
		}

		if (fname != null) // a seal configuration file was specified
		{
			if (!fname.exists())
				throw new ParseException("Configuration file " + fname + " doesn't exist");
			if (!fname.canRead())
				throw new ParseException("Can't read configuration file " + fname);
			// at this point it should be all good.
		}
		else // none specified.  Try the default
		{
			// presume that if it exists the user intends to use it
			if (DefaultConfigFile.exists())
			{
				if (DefaultConfigFile.canRead())
					fname = DefaultConfigFile;
				else
				{
					// The file exists but it can't be read.  Warn the user.
					LogFactory.getLog(SealToolParser.class).warn("Seal configuration file " + DefaultConfigFile + " isn't readable");
					// leave fname as null so no configuration file will be used
				}
			}
		}

		return fname;
	}

	/**
	 * Set properties useful for the whole Seal suite.
	 */
	protected void setDefaultProperties(Configuration conf)
	{
		conf.set("mapred.compress.map.output", "true");
	}

	public void parse(Configuration conf, String[] args) throws IOException
	{
		try
		{
			parseOptions(conf, args);
		}
		catch( ParseException e )
		{
			defaultUsageError(e.getMessage()); // doesn't return
		}
	}

	/**
	 * Parses command line.
	 *
	 * Override this method to implement additional command line options,
	 * but do make sure you call this method to parse the default options.
	 */
	protected CommandLine parseOptions(Configuration conf, String[] args)
	  throws ParseException, IOException
	{
		myconf = conf;

		setDefaultProperties(conf);

		// load settings from configuration file
		// first, parse the command line (in getRcFile) looking for an option overriding the default seal configuration file
		File configFile = getRcFile(args);
		if (configFile != null)
			loadConfig(conf, configFile);

		// now parse the entire command line using the default hadoop parser.  Now
		// the user can override properties specified in the config file with properties
		// specified on the command line.
		CommandLine line = new GenericOptionsParser(conf, options, args).getCommandLine();
		if (line == null)
			throw new ParseException("Error parsing command line"); // getCommandLine returns an null if there was a parsing error

		////////////////////// input/output formats //////////////////////
		// set the configuration property.  Then, we'll check the property
		// to ensure it has a valid value, regardless of whether we just set it,
		// so that the check will also be valid if the property is set directly.
		if (line.hasOption(opt_inputFormat.getOpt()))
			myconf.set(INPUT_FORMAT_CONF, line.getOptionValue(opt_inputFormat.getOpt()));

		validateIOFormat(INPUT_FORMAT_CONF, acceptedInputFormats);

		if (line.hasOption(opt_outputFormat.getOpt()))
			myconf.set(OUTPUT_FORMAT_CONF, line.getOptionValue(opt_outputFormat.getOpt()));

		validateIOFormat(OUTPUT_FORMAT_CONF, acceptedOutputFormats);

		/////////////////////// output compression /////////////////////
		if (line.hasOption(opt_compressOutput.getOpt()))
		{
			myconf.setBoolean("mapred.output.compress", true);
			String codec = line.getOptionValue(opt_compressOutput.getOpt());
			if (codec != null)
			{
				String codecClass = "org.apache.hadoop.io.compress.GzipCodec"; // default
				if ("auto".equalsIgnoreCase(codec) || "gzip".equalsIgnoreCase(codec))
				{
					// pass.  Already set
				}
				else if ("bzip2".equalsIgnoreCase(codec))
					codecClass = "org.apache.hadoop.io.compress.BZip2Codec";
				else if ("snappy".equalsIgnoreCase(codec))
					codecClass = "org.apache.hadoop.io.compress.SnappyCodec";
				else
				{
					throw new ParseException("Unknown codec " + codec + ". Valid values are gzip, bzip2, snappy and auto.\n" +
					   "If you want to use an unsupported codec pass 'auto' and set the property mapred.output.compression.codec directly");
				}

				myconf.set("mapred.output.compression.codec", codecClass);
			}
		}


		////////////////////// number of reducers //////////////////////
		if (line.hasOption(opt_nReduceTasks.getOpt()))
		{
			String rString = line.getOptionValue(opt_nReduceTasks.getOpt());
			try
			{
				int r = Integer.parseInt(rString);
				if (r >= minReduceTasks)
					nReduceTasks = r;
				else
					throw new ParseException("Number of reducers must be greater than or equal to " + minReduceTasks + " (got " + rString + ")");
			}
			catch (NumberFormatException e)
			{
				throw new ParseException("Invalid number of reduce tasks '" + rString + "'");
			}
		}

		////////////////////// positional arguments //////////////////////
		String[] otherArgs = line.getArgs();
		if (otherArgs.length < 2) // require at least two:  one input and one output
			throw new ParseException("You must provide input and output paths");
		else
		{
			//
			FileSystem fs;
			for (int i = 0; i < otherArgs.length - 1; ++i) {
				Path p = new Path(otherArgs[i]);
				fs = p.getFileSystem(conf);
				p = p.makeQualified(fs);
				FileStatus[] files = fs.globStatus(p);
				if (files != null && files.length > 0)
				{
					for (FileStatus status: files)
						inputs.add(status.getPath());
				}
				else
					throw new ParseException("Input path " + p.toString() + " doesn't exist");
			}
			// now the last one, should be the output path
			outputDir = new Path(otherArgs[otherArgs.length - 1]);
			fs = outputDir.getFileSystem(conf);
			outputDir = outputDir.makeQualified(fs);
			if (fs.exists(outputDir))
				throw new ParseException("Output path " + outputDir.toString() + " already exists.  Won't overwrite");
		}

		return line;
	}

	/**
	 * Set accepted input format names.
	 *
	 * @names if null, turns off checking of input format name.  Otherwise,
	 *        specifying as a parameter to --input-format any name that is
	 *        not in the list will generate an error.
	 */
	public void setAcceptedInputFormats(String[] names)
	{
		acceptedInputFormats = names;

		String inputFormatHelp = INPUT_FORMAT_DESC;
		if (names != null)
			inputFormatHelp += " (" + joinStrings(names, ",") + ")";

		opt_inputFormat.setDescription(inputFormatHelp);
	}

	/**
	 * Set accepted output format names.
	 *
	 * @names if null, turns off checking of output format name.  Otherwise,
	 *        specifying as a parameter to --output-format any name that is
	 *        not in the list will generate an error.
	 */
	public void setAcceptedOutputFormats(String[] names)
	{
		acceptedOutputFormats = names;

		String help = OUTPUT_FORMAT_DESC;
		if (names != null)
			help += "(" + joinStrings(names, ",") + ")";

		opt_outputFormat.setDescription(help);
	}

	/**
	 * Return the input format specified, if any.
	 *
	 * @param defaultName Return this value if the configuration isn't set.
	 */
	public String getInputFormatName(String defaultName)
	{
		return myconf.get(INPUT_FORMAT_CONF, defaultName);
	}

	/**
	 * Return the input format specified, if any.
	 */
	public String getInputFormatName()
	{
		return myconf.get(INPUT_FORMAT_CONF);
	}

	/**
	 * Return the output format specified, if any.
	 *
	 * @param defaultName Return this value if the configuration isn't set.
	 */
	public String getOutputFormatName(String defaultName)
	{
		return myconf.get(OUTPUT_FORMAT_CONF, defaultName);
	}

	/**
	 * Return the output format specified, if any.
	 */
	public String getOutputFormatName()
	{
		return myconf.get(OUTPUT_FORMAT_CONF);
	}

	/**
	 * Get total number of reduce tasks to run.
	 * This option parser must have already parsed the command line.
	 */
	public int getNReduceTasks() throws java.io.IOException
 	{
		if (myconf == null)
			throw new IllegalStateException("getNReduceTasks called before parsing the command line.");

		if (nReduceTasksPerNode < 0)
			throw new IllegalArgumentException("Invalid number of default reduce tasks per node: " + nReduceTasksPerNode);

		if (nReduceTasks == null)
		{
			// Calculate and cache value
			// To calculate, we get the number of tasktrackers available and multiply by
			// nReduceTasksPerNode.  It can happen that a cluster doesn't yet have any
			// task trackers, so we use a lower bound of 1.
			nReduceTasks = Math.max(ClusterUtils.getNumberTaskTrackers(myconf), 1) * nReduceTasksPerNode;
			return nReduceTasks;
		}
		else
			return nReduceTasks;
 	}

	public void setNReduceTasksPerNode(int value)
	{
		if (value < 0)
			throw new IllegalArgumentException("number of reduce tasks per node must be >= 0 (got " + value + ")");
		nReduceTasksPerNode = value;
		nReduceTasks = null; // reset cached value
	}

	/**
	 * Return the specified output path.
	 */
	public Path getOutputPath()
	{
		return outputDir;
	}

	public List<Path> getInputPaths()
	{
		ArrayList<Path> retval = new ArrayList<Path>(getNumInputPaths());
		for (Path p: inputs)
			retval.add(p);
		return retval;
	}

	public int getNumInputPaths()
	{
		return inputs.size();
	}

	public void defaultUsageError()
	{
		defaultUsageError(null);
	}

	/**
	 * Prints help and exits with code 3.
	 */
	public void defaultUsageError(String msg)
	{
		System.err.print("Usage error");
		if (msg != null)
			System.err.println(":  " + msg);
		System.err.print("\n");
		// XXX: redirect System.out to System.err since the simple version of
		// HelpFormatter.printHelp prints to System.out, and we're on a way to
		// a fatal exit.
		System.setOut(System.err);
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("" + toolName + " [options] <in>+ <out>", options);
		System.exit(3);
	}

	protected static String joinStrings(String[] strings, String joinString)
	{
		if (strings.length == 0)
			return "";
		else if (strings.length == 1)
			return strings[0];
		else
		{
			StringBuilder builder = new StringBuilder();
			builder.append(strings[0]);
			for (int i = 1; i < strings.length; ++i)
			{
				builder.append(joinString);
				builder.append(strings[i]);
			}

			return builder.toString();
		}
	}

	protected void validateIOFormat(String ioProperty, String[] acceptedFormats) throws ParseException
	{
		String selectedFormat = myconf.get(ioProperty);

		if (acceptedFormats != null && selectedFormat != null)
		{
			for (int i = 0; i < acceptedFormats.length; ++i)
				if (acceptedFormats[i].equals(selectedFormat))
					return;
			throw new ParseException("Incompatible file format selected. " + ioProperty + " is " + selectedFormat + " but acceptable values are " + joinStrings(acceptedFormats, ", "));
		}
	}
}
