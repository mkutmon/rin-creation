package cytargetlinker.conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import org.bridgedb.DataSource;
import org.bridgedb.IDMapperException;
import org.bridgedb.IDMapperStack;
import org.bridgedb.Xref;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.utils.ArgsParser;
import cytargetlinker.conversion.utils.CommonAttributes;
import cytargetlinker.conversion.utils.Utils;
import cytargetlinker.conversion.utils.ArgsParser.AFilesAttributes;
import cytargetlinker.conversion.utils.ArgsParser.AFilesIn;
import cytargetlinker.conversion.utils.ArgsParser.AFilesOut;
import cytargetlinker.conversion.utils.ArgsParser.AHelp;
import cytargetlinker.conversion.utils.ArgsParser.GraphBuilder;

/**
 * download the TF-gene networks from
 * http://encodenets.gersteinlab.org/enets2.Proximal_filtered.txt
 * http://encodenets.gersteinlab.org/enets3.Distal.txt
 * 
 * you can find the annotation file (GENCODE) here:
 * http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeGencodeV12/wgEncodeGencodeAttrsV12.tab.gz
 * 
 * the family info can be found in a pdf (supplementary data)
 * for the annotation create a file with three columns: TF / TF Class / TF Family
 * 
 * @author martina
 *
 */
public class EncodeNetworks {

	private final static Logger log = Logger.getLogger(MiRecords.class.getName());
	private static Args pargs;
	private interface Args extends AHelp, AFilesIn, AFilesOut, AFilesAttributes {}
	
	/**
	 * USAGE: java -cp conversion.jar cytargetlinker.conversion.EncodeNetworks 
	 * ARGUMENTS: 
	 * -i = input file
	 * [-o] = output file
	 * [--bridgeDbFile] = BridgeDb mapping files for source and target nodes
	 * [-a] = Annotation file (GENCODE attributes)
	 * [--organism] = "Homo sapiens"
	 * [-f] = TF family file
	 */
	public static void main(String argv[]) throws Exception {
		pargs = ArgsParser.parse(argv, Args.class);
		if(pargs.getAnnotationFile() != null) {
			EncodeNetworks converter = new EncodeNetworks();
			converter.startConversion();
		} else {
			log.severe("Please add annotation file to map from Ensembl transcript to Ensembl gene id.");
		}
	}
	
	private IDMapperStack idMapper;
	private Graph graph;
	private Map<String, List<String>> edges;
	
	private List<String> foundConnections;
	private Map<String, String> annotationMap;
	private static Map<String, List<String>> annotationErrors;
	private Map<String, String> familyInfo;
	private List<String> tfs;
	private List<String> targets;
	
	private int countEdge = 0;
	private int countGenes = 0;
	private int countTFs = 0;
	
	private boolean mapping = false;

	public EncodeNetworks() throws Exception {
		edges = new HashMap<String, List<String>>();
		annotationMap = new HashMap<String, String>();
		foundConnections = new ArrayList<String>();
		annotationErrors = new HashMap<String, List<String>>();
		familyInfo = new HashMap<String, String>();
		tfs = new ArrayList<String>();
		targets = new ArrayList<String>();
	}

	private void startConversion() throws Exception {
		File in = pargs.getInput();
		if(in != null) {
			if(pargs.getAnnotationFile() != null) {
				Utils.setUpLogger(log, getLogFile(), false);

				if(pargs.isBridgeDbFiles()) {
					idMapper = Utils.initIDMapper(pargs.getBridgeDbFiles(), false);
					if(idMapper != null) {
						mapping = true;
						log.info("identifier mapping on");
					}
				}
				
				log.info("conversion of encode networks started ...\n");
	
				ArgsParser.convertAndWrite(pargs, pargs, new GraphBuilder() {
					public Graph buildGraph(File in) throws Exception {
						return convert(in);
					}
				});
					
				log.info("Encode networks file created.\n\n");
			} else {
				log.severe("Please specify annotation file (GENCODE).");
			}
		} else {
			log.severe("Please specify input file!");
		}
	}
	
	
	
	private Graph convert(File in) {
		graph = new Graph();
		try {
			readAnnotations();
			if(pargs.isTFFamilyFile()) {
				readFamilyInfo(pargs.getTFFamilyFile());
			}
			
			setNetworkAttributes(in);
			
			try {
				BufferedReader br = new BufferedReader(new FileReader(in));
				String line;
				while((line = br.readLine()) != null) {
					String [] buffer = line.split("\t");
					if(buffer.length != 3) {
						buffer = line.split(" ");
					}
					String source = buffer[0];
					String target = buffer[2];
						
					if(!tfs.contains(source)) {
						tfs.add(source);
						if(targets.contains(source)) {
							targets.remove(source);
						}
					}
					if(!targets.contains(target)) {
						if(!tfs.contains(target)) {
							targets.add(target);
						}
					}
				}
				br.close();
			} catch (IOException e) {
				log.warning("Could not read input file " + in.getAbsolutePath());
			}
			try {		
				BufferedReader br = new BufferedReader(new FileReader(in));
				String line;
				int mapped = 0;
				int notMapped = 0;
				while((line = br.readLine()) != null) {
					String [] buffer = line.split("\t");
					if(buffer.length != 3) {
						buffer = line.split(" ");
					}
					String source = buffer[0];
					String intType = buffer[1];
					String target = buffer[2];
						
					if(annotationMap.containsKey(source) && annotationMap.containsKey(target)) {
						createNode(annotationMap.get(source), source, "transcriptionFactor");
						if(tfs.contains(target)) {
							createNode(annotationMap.get(target), target, "transcriptionFactor");
						} else {
							createNode(annotationMap.get(target), target, "gene");
						}
							
						createEdge(annotationMap.get(source), annotationMap.get(target), intType);
						mapped++;
					} else {
						if(annotationMap.containsKey(source)) log.info("\t target not mapped in annotation file: " + target);
						else log.info("\t source not mapped in annotation file: " + source);
						notMapped++;
					}
				}
					
				br.close();
				log.info(mapped + " interactions were mapped.\n" + notMapped + " interactions were not mapped.");
			} catch (IOException e) {
				log.warning("Could not read input file " + in.getAbsolutePath());
			}
			log.info(foundConnections.size() + " interactions have been found.\n" + countGenes + " gene nodes.\n" + countTFs + " transcription factor nodes.\n");
		} catch (IOException e) {
			log.severe("Could not parse annotation file.");
		}
		return graph;
	}
	
	private void createEdge(String source, String target, String intType) {
		if(edges.containsKey(source)) {
			if(!edges.get(source).contains(target)) {
				Edge e = graph.addEdge("" + countEdge, graph.getNode(source),graph.getNode(target));
				e.setAttribute("datasource", "ENCODE network (" + pargs.getDescription() + ")");
				e.setAttribute("interactionType", intType);
				
				edges.get(source).add(target);
				foundConnections.add(source + "\t" + target);
				countEdge++;
			}
		} else {
			Edge e = graph.addEdge("" + countEdge, graph.getNode(source),graph.getNode(target));
			e.setAttribute("datasource", "ENCODE network (" + pargs.getDescription() + ")");
			e.setAttribute("interactionType", intType);
			
			List<String> list = new ArrayList<String>();
			list.add(target);
			foundConnections.add(source + "\t" + target);
			edges.put(source, list);
			countEdge++;
		}
	}

	private void createNode(String id, String name, String type) {
		
		if(graph.getNode(id) == null) {
			String entrez = "";
			String identifiers = "[" + id;
			
			if(mapping) {
				Xref xrefIn = new Xref(id, DataSource.getBySystemCode("En"));
				Set<Xref> e;
				try {
					e = idMapper.mapID(xrefIn, DataSource.getBySystemCode("L"));
					if(!e.isEmpty()) {
						entrez = e.iterator().next().getId();
						identifiers = identifiers + "," + entrez;
					}
				} catch (IDMapperException ex) {
					// could not be mapped
				}
				
			}
			identifiers = identifiers + "]";
			
			Node node = graph.addNode(id);
			node.appendAttribute("identifiers", identifiers);
			node.appendAttribute("label", name);
			node.appendAttribute("name", name);
			node.appendAttribute("biologicalType", type);
			node.appendAttribute("ensemblID", id);
			node.appendAttribute("entrezGeneID", entrez);
			
			if(familyInfo.containsKey(name)) {
				String [] buffer = familyInfo.get(name).split(";");
				node.appendAttribute("TF-Class", buffer[0]);
				node.appendAttribute("TF-Family", buffer[1]);
			}
			
			if(type.equals("gene")) {
				countGenes++;
			} else {
				countTFs++;
			}
		} else {
			Node node = graph.getNode(id);
			if(!node.getAttribute("biologicalType").equals(type)) {
				if(type.equals("transcriptionFactor")) {
					node.setAttribute("biologicalType", type);
				}
			}
		}
	}
	
	private void setNetworkAttributes(File in) {
		graph.setTitle("ENCODE network (" + pargs.getDescription() + ")");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), "ENCODE network (" + pargs.getDescription() + ")");
	
		graph.setAttribute(CommonAttributes.TYPE.getName(), "TF-gene interactions");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), in.getName());
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), "Ensembl");
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "Ensembl");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "transcriptionFactor");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "gene");
	}

	private void readAnnotations() throws IOException {
		File file = new File(pargs.getAnnotationFile());
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		// read header line
		reader.readLine();
		
		String line;
		while((line = reader.readLine()) != null) {
			String [] str = line.split("\t");
			String name = str[1];
			String id = str[0];
			
			if(id.contains(".")) {
				int index = id.indexOf(".");
				id = id.substring(0,index);
			}
			
			if(annotationMap.containsKey(name)) {
				if(!annotationMap.get(name).equals(id)) {
					if(annotationErrors.containsKey(name)) {
						annotationErrors.get(name).add(id);
					} else {
						List<String> l = new ArrayList<String>();
						l.add(id);
						annotationErrors.put(name, l);
					}
				}
			} else {
				annotationMap.put(name, id);
			}
		}
		
		// remove gene names with more than one id
		for(String key : annotationErrors.keySet()) {
			annotationMap.remove(key);
		}
		
		reader.close();
	}
	
	private void readFamilyInfo(File file) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(file));
		reader.readLine();
		String line;
		while((line = reader.readLine()) != null) {
			String [] str = line.split("\t");
			String name = str[0];
			String tfClass = str[1];
			String tfFamily = str[2];
			
			if(!familyInfo.containsKey(name)) {
				familyInfo.put(name, tfClass + ";" + tfFamily);
			}
		}
	}
	
	private File getLogFile() {
		if(pargs.isLogFile()) {
			return pargs.getLogFile();
		} else {
			if(pargs.isInput()) {
				return new File(pargs.getInput() + ".log");
			} else {
				return null;
			}
		}
	}
}