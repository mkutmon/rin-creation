package cytargetlinker.conversion.utils;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.util.logging.Logger;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;
import uk.co.flamingpenguin.jewel.cli.Option;
import cytargetlinker.conversion.graph.GmlWriter;
import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.XGMMLWriter;

/**
 * Commonly used command line arguments to be parsed.
 * @author Thomas
 */
public class ArgsParser {
	private final static Logger log = Logger.getLogger(ArgsParser.class.getName());
	
	public static <A> A parse(String[] args, Class<? extends A> c) throws ArgumentValidationException {
		return CliFactory.parseArguments(c, args);
	}

	public interface AHelp {
		@Option(helpRequest = true, shortName = "h")
		boolean getHelp();
	}
	
	public interface AFilesIn {
		@Option(shortName = "i", description = "The path to the input file.")
		public File getInput();
		public boolean isInput();
	}
	
	public interface AFilesOut {
		@Option(shortName = "o", description = "The output GML file or directory to write the network(s) to.")
		public File getOutput();
		public boolean isOutput();
	}

	public interface AFilesAttributes {
		@Option(shortName = "n", description = "The name of the network.")
		public String getName();
		public boolean isName();
		
		@Option(shortName = "a", description = "Annotation file")
		public String getAnnotationFile();
		public boolean isAnnotationFile();
		
		@Option(longName = "bridgeDbFile", description = "BridgeDb mapping files")
		public List<File> getBridgeDbFiles();
		public boolean isBridgeDbFiles();
		
		@Option(longName = "organism", description = "BridgeDb mapping file")
		public String getOrganism();
		public boolean isOrganism();
		
		@Option(shortName = "s", description = "TFe sources (auto or user)")
		public String getSource();
		public boolean isSource();
		
		@Option(shortName = "l", description = "The path to the log file.")
		public File getLogFile();
		public boolean isLogFile();
		
		@Option(shortName = "d", description = "Description.")
		public String getDescription();
		public boolean isDescription();
		
		@Option(shortName = "f", description = "TF Family file")
		public File getTFFamilyFile();
		public boolean isTFFamilyFile();
		
		@Option(shortName = "v", description = "database version")
		public String getDatabaseVersion();
		public boolean isDatabaseVersion();
	}
	
	public interface AIDMapper {
		@Option(description = "Bridgedb connection strings to idmappers to use to translate between xrefs of different types.")
		public List<String> getIdm();
		public boolean isIdm();
		
		@Option(defaultValue = { "L", "Ce" }, description = "The datasource(s) to translate all xrefs to (use system code).")
		public List<String> getDs();
	}
	
	public interface GraphBuilder {
		public Graph buildGraph(File in) throws Exception;
	}
	
	private interface GraphWriter {
		public void write(Graph g, PrintWriter out) throws Exception;
	}
	
	private static class XGMML implements GraphWriter {
		public void write(Graph g, PrintWriter out) throws Exception { XGMMLWriter.write(g, out); }
	}
	
	private static class GML implements GraphWriter {
		public void write(Graph g, PrintWriter out) throws Exception { GmlWriter.write(g, out); }
	}
	
	/**
	 * writes xgmml file
	 * @param fi
	 * @param fo
	 * @param gb
	 * @throws Exception
	 */
	public static void convertAndWrite(AFilesIn fi, AFilesOut fo, GraphBuilder gb) throws Exception {
		File input = fi == null ? null : fi.getInput();
		File output = fo.isOutput() ? fo.getOutput() : new File(fi.getInput().getAbsolutePath() + ".xgmml");
		
		log.info("Converting " + input + " to " + output + "\n");
		
		Graph g = gb.buildGraph(input);
		g.setDirected(true);
		
		GraphWriter writer = new XGMML();
		if(output.getName().endsWith(".gml")) {
			writer = new GML();
		}
		PrintWriter po = new PrintWriter(output);
		writer.write(g, po);
		po.close();
	}
}