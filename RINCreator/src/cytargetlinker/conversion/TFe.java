package cytargetlinker.conversion;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.logging.Logger;

import org.apache.commons.io.IOUtils;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.utils.ArgsParser;
import cytargetlinker.conversion.utils.CommonAttributes;
import cytargetlinker.conversion.utils.Utils;
import cytargetlinker.conversion.utils.ArgsParser.AFilesAttributes;
import cytargetlinker.conversion.utils.ArgsParser.AFilesOut;
import cytargetlinker.conversion.utils.ArgsParser.AHelp;
import cytargetlinker.conversion.utils.ArgsParser.GraphBuilder;

/**
 * Import transcription factor targets from the TFe wiki.
 * @author thomas
 */
public class TFe {
	private final static Logger log = Logger.getLogger(TFe.class.getName());

	private static Args arguments;
	private interface Args extends AHelp, AFilesOut, AFilesAttributes { }
	
	private final String TFE_URL = "http://www.cisreg.ca/cgi-bin/tfe/api.pl?";
	
	private int numInteractions = 0;
	private Graph graph;
	
	public TFe() {
		
	}
	
	public void startConversion() throws Exception {
		ArgsParser.convertAndWrite(null, arguments, new GraphBuilder() {
			public Graph buildGraph(File in) throws Exception {
				return importTFe(null);
			}
		});
	}
	
	
	public Graph importTFe(String species) throws MalformedURLException, IOException {
		if(getLogFile() != null) {
			Utils.setUpLogger(log, getLogFile(), false);
			String type;
			if(!arguments.isSource()) {
				type = "both";
				createGraph(species, type);
			} else {
				if(arguments.getSource().equals("user") || arguments.getSource().equals("auto")) {
					createGraph(species, arguments.getSource());
				} else {
					System.out.println("Please specify a correct TFe source type [user or auto]\t" + arguments.getSource());
				}
			}
		} else {
			System.out.println("Please specify log file!");
		}

		return graph;
	}
	
	private void createGraph(String species, String type) throws MalformedURLException, IOException {
		graph = new Graph();
		graph.setTitle("TFe");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), "TFe");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), TFE_URL);
		graph.setAttribute(CommonAttributes.TYPE.getName(), "TF targets");
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), "Entrez Gene");
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "Entrez Gene");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "gene");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "gene");
		
		//First get all TF ids
		String[] tfIds = readURL(TFE_URL + "code=all-tfids").split("\n");
		for(String id : tfIds) {
			log.info("Processing " + id);
			if(species != null) {
				String tfSpecies = readURL(TFE_URL + "tfid=" + id + "&code=species").trim();
				log.info(tfSpecies);
				if(!species.equals(tfSpecies)) continue;
			}
	
			log.info("Querying info for " + id);
			String tfId = readURL(TFE_URL + "tfid=" + id + "&code=entrez-gene-id").trim();
			
			//Get targets
			String[] lines = readURL(TFE_URL + "tfid=" + id + "&code=targets").split("\n");
			if(lines.length > 0) {
				boolean created = false;
				Node nSrc = null;
				for(String l : lines) {
					if("".equals(l)) continue;
					String[] cols = l.split("\t");
					for(int i = 0; i < cols.length; i++) cols[i] = cols[i].trim();
					
					String targetId = cols[0];

					if(type.equals("both") || cols[5].equals(type)) {
						if(!created) {
							nSrc = createSourceNode(tfId, id, graph);
							created = true;
						}
						Node nTgt = graph.addNode(targetId);
						nTgt.setAttribute("Label", cols[1]);
						nTgt.setAttribute("EntrezId", targetId);
						nTgt.appendAttribute("TFeActingComplex", cols[2]);
						
						Edge edge = graph.addEdge(nSrc + "target" + nTgt, nSrc, nTgt);
						nTgt.setAttribute("TFeActingComplex", cols[2]);
						edge.setAttribute("TFeEffect", cols[3]);
						edge.setAttribute("PMID", cols[4]);
						edge.setAttribute("TFeSource", cols[5]);
						edge.setAttribute("Interaction", "TF target");
						edge.setAttribute("Directed", "true");
						numInteractions++;
					} 
				}
			}
		}
		
		System.out.println("Number of interactions (type = " + type + "): " + numInteractions);
	}

	private Node createSourceNode(String tfId, String id, Graph graph) throws MalformedURLException, IOException {
		Node nSrc = graph.addNode(tfId);
		nSrc.setAttribute("Label", readURL(TFE_URL + "tfid=" + id + "&code=symbol").trim());
		nSrc.setAttribute("EntrezId", tfId);
		nSrc.setAttribute("Species", readURL(TFE_URL + "tfid=" + id + "&code=species").trim());
		return nSrc;
	}
	
	private String readURL(String url) throws MalformedURLException, IOException {
		return IOUtils.toString(new URL(url).openStream());
	}
	
	public static void main(String[] args) {
		try {
			arguments = ArgsParser.parse(args, Args.class);
			TFe converter = new TFe();
			converter.startConversion();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private File getLogFile() {
		if(arguments.isLogFile()) {
			return arguments.getLogFile();
		} else {
			return null;
		}
	}
}
