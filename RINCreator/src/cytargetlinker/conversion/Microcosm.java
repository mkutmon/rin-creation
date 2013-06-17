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

import cytargetlinker.conversion.data.GeneNode;
import cytargetlinker.conversion.data.MTI;
import cytargetlinker.conversion.data.MiRNANode;
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
 * Converts miRNA target text files from Microcosm to a XGMML or GML network.
 * 
 * The input files for this script can be downloaded here:
 * 
 * @author Thomas, Pooja, Naime, Tina
 */
public class Microcosm {
	
	private final static Logger log = Logger.getLogger(MiRecords.class.getName());
	private static Args pargs;
	private interface Args extends AHelp, AFilesIn, AFilesOut, AFilesAttributes {}
	
	/**
	 * MAIN METHOD
	 * USAGE: java -cp conversion.jar cytargetlinker.conversion.Microcosm 
	 * ARGUMENTS: 
	 * -i = input file
	 * -o = output file
	 * --bridgeDbFile = BridgeDb mapping files for source and target nodes
	 * -a = Annotation file (csv containing Transcript/Gene mapping - separated with a comma)
	 * --organism "mouse" or "human" 
	 */
	public static void main(String argv[]) throws Exception {
		pargs = ArgsParser.parse(argv, Args.class);
		if(pargs.getAnnotationFile() != null) {
			Microcosm converter = new Microcosm();
			converter.startConversion();
		} else {
			log.severe("Please add annotation file to map from Ensembl transcript to Ensembl gene id.");
		}
	}
	
	private IDMapperStack gdb;
	private Graph graph;
	private Map<String, Integer> index;
	private Map<String, List<String>> edges;
	
	private List<String> genesNotFound;
	private List<String> foundConnections;
	private Map<String, String> annotationMap;
	
	private int countEdge = 0;
	private boolean mapping = true;
	private int countGenes = 0;
	private int countMiRNAs = 0;
	
	private Map<String, MiRNANode> miRNAs;
	private Map<String, GeneNode> genes;
	private Map<String, MTI> interactions;
	
	private List<String> nodesNotFound;
	private List<String> interactionIgnored;
	
	public Microcosm() throws Exception {
		edges = new HashMap<String, List<String>>();
		index = new HashMap<String, Integer>();
		annotationMap = new HashMap<String, String>();
		foundConnections = new ArrayList<String>();
		genesNotFound = new ArrayList<String>();
		
		miRNAs = new HashMap<String, MiRNANode>();
		genes = new HashMap<String, GeneNode>();
		interactions = new HashMap<String, MTI>();
		nodesNotFound = new ArrayList<String>();
		interactionIgnored = new ArrayList<String>();
	}
	
	private void startConversion() throws Exception {
		File in = pargs.getInput();
		if(in != null) {
			Utils.setUpLogger(log, getLogFile(), false);
			if(pargs.isBridgeDbFiles()) gdb = Utils.initIDMapper(pargs.getBridgeDbFiles(), false);
			if(gdb == null) {
				mapping = false;
				log.info("no identifier mapping");
			}
			log.info("conversion of microcosm file started ...\n");

			ArgsParser.convertAndWrite(pargs, pargs, new GraphBuilder() {
				public Graph buildGraph(File in) throws Exception {
					return importMicrocosm(in);
				}
			});
				
			log.info("MicroCosm network file created.\n\n");
		} else {
			System.out.println("Please specify input file!");
		}
	}
	
	protected Graph importMicrocosm(File in) throws IOException {
		graph = new Graph();
		
		readAnnotations();
		setNetworkAttributes(in);
		parseMicrocosm(in);
		
		log.info(foundConnections.size() + " interactions have been found.\n" + genesNotFound.size() + " transcripts were not mapped to genes.\n" + countGenes + " gene nodes.\n" + countMiRNAs + " miRNA nodes.\n");
		cleanUp();
		return graph;
	}
	
	private void cleanUp() {
		edges.clear();
		index.clear();
		annotationMap.clear();
		foundConnections.clear();
		genesNotFound.clear();
		
		genes.clear();
		miRNAs.clear();
		interactions.clear();
		nodesNotFound.clear();
		interactionIgnored.clear();
	}
	
	private List<String> ids;
	
	private void parseMicrocosm(File in) {
		ids = new ArrayList<String>();

		try {
			BufferedReader br = new BufferedReader(new FileReader(in));
			String[] header = readHeader(br);

			if (header != null) {
				String line = null;

				List<String[]> rows = new ArrayList<String[]>();
				while ((line = br.readLine()) != null) {
					if (!line.startsWith("#") && !line.equals("")) {
						rows.add(line.split("\t", header.length));
					}
				}

				for (String[] row : rows) {
					String transcriptId = row[index.get("TRANSCRIPT_ID")];
					String geneName = row[index.get("EXTERNAL_NAME")];
					String mirna = row[index.get("SEQ")];
					String score = row[index.get("SCORE")];
					String pValue = row[index.get("PVALUE_OG")];

					if (annotationMap.containsKey(transcriptId)) {
						if (!ids.contains(annotationMap.get(transcriptId))) {
							createGeneNode(annotationMap.get(transcriptId),
									geneName);
							ids.add(annotationMap.get(transcriptId));
						}
						if (!ids.contains(mirna)) {
							createMiRNANode(mirna);
							ids.add(mirna);
						}
						addEdge(annotationMap.get(transcriptId), mirna, score,
								pValue);
					} else {
						if (!genesNotFound.contains(transcriptId))
							genesNotFound.add(transcriptId);
					}
				}
			}
		} catch (IOException e) {
			log.warning("Could not read input file " + in.getAbsolutePath());
		}
	}
	
	private void addEdge(String gene, String mirna, String score, String pValue) {
		if(edges.containsKey(gene)) {
			if(!edges.get(gene).contains(mirna)) {
				Edge e = graph.addEdge("" + countEdge, graph.getNode(mirna),graph.getNode(gene));
				e.setAttribute("datasource", "Microcosm Targets version 5");
				e.setAttribute("interactionType", "predicted MTI");
				try {
					Double d = Double.parseDouble(score);
					e.setAttribute("score", d.toString());
					Double pvalue = Double.parseDouble(pValue);
					e.setAttribute("pvalue", pvalue.toString());
				} catch(NumberFormatException ex) {
					e.setAttribute("score", "");
					e.setAttribute("pvalue", "");	
				}
				
				edges.get(gene).add(mirna);
				foundConnections.add(gene + "\t" + mirna);
				countEdge++;
			}
		} else {
			Edge e = graph.addEdge("" + countEdge, graph.getNode(mirna),graph.getNode(gene));
			e.setAttribute("datasource", "Microcosm Targets version 5");
			e.setAttribute("interactionType", "predicted MTI");
			try {
				Double d = Double.parseDouble(score);
				e.setAttribute("score", d.toString());
				Double pvalue = Double.parseDouble(pValue);
				e.setAttribute("pvalue", pvalue.toString());
			} catch(NumberFormatException ex) {
				e.setAttribute("score", "");
				e.setAttribute("pvalue", "");
			}
			List<String> list = new ArrayList<String>();
			list.add(mirna);
			foundConnections.add(gene + "\t" + mirna);
			edges.put(gene, list);
			countEdge++;
		}
	}
	
	private String[] readHeader(BufferedReader br) throws IOException {
		String[] header = null;
		String line = null;
		while((line = br.readLine()) != null) {
			if(line.startsWith("##GROUP")) {
				header = line.split("\t");
				index = new HashMap<String, Integer>();
				
				for(int i = 0; i < header.length; i++) index.put(header[i], i);
				break;
			}
		}
		return header;
	}
	
	private void setNetworkAttributes(File input) {
		graph.setTitle("Microcosm " + pargs.getOrganism());
		graph.setAttribute(CommonAttributes.DATABASE.getName(), "Microcosm (" + pargs.getOrganism() + ")");
	
		graph.setAttribute(CommonAttributes.TYPE.getName(), "predicted miRNA targets");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), input.getName());
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), "Ensembl");
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "miRBase");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "mirna");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "gene");
	}

	private void readAnnotations() {
		File file = new File(pargs.getAnnotationFile());
		if(file.exists()) {
			try {
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String line;
				while((line = reader.readLine()) != null) {
					String [] buffer = line.split(",");
					if(!annotationMap.containsKey(buffer[1])) {
						annotationMap.put(buffer[1], buffer[0]);
					}
				}
				reader.close();
			} catch (IOException e) {
				log.severe("Could not read annotation file.\n");
				System.exit(0);
			}
		} else {
			log.severe("Annotation file does not exist.\n");
			System.exit(0);
		}
	}

	private void createGeneNode(String geneId, String geneName) {		
		String type = "gene";
		String organism = pargs.getOrganism();
		
		if(graph.getNode(geneId) == null) {
			String entrez = "";
			String identifiers = "[" + geneId;
			
			if(mapping) {
				Xref xrefIn = new Xref(geneId, DataSource.getBySystemCode("En"));
				Set<Xref> e;
				try {
					e = gdb.mapID(xrefIn, DataSource.getBySystemCode("L"));
					if(!e.isEmpty()) {
						entrez = e.iterator().next().getId();
						identifiers = identifiers + "," + entrez;
					}
				} catch (IDMapperException ex) {
					// could not be mapped
				}
				
			}
			identifiers = identifiers + "]";
			
			Node node = graph.addNode(geneId);
			node.appendAttribute("identifiers", identifiers);
			node.appendAttribute("label", geneName);
			node.appendAttribute("biologicalType", type);
			node.appendAttribute("organism", organism);
			node.appendAttribute("ensemblID", geneId);
			node.appendAttribute("entrezGeneID", entrez);
			countGenes++;
		}		
	}

	private void createMiRNANode(String miRNA) {
		String type = "microRNA";
		String organism = pargs.getOrganism();
		
		if(graph.getNode(miRNA) == null) {
			String mimat = "";
			String identifiers = "[" + miRNA;
			Xref xrefIn = new Xref(miRNA, DataSource.getBySystemCode("Mb"));
			
			if(mapping) {
				try {
					Set<Xref> result = gdb.mapID(xrefIn, DataSource.getBySystemCode("Mb"));
					Set<Xref> result2 = gdb.mapID(xrefIn, DataSource.getBySystemCode("Mbm"));
					
					List<String> list = new ArrayList<String>();
					list.add(miRNA);
					for(Xref x : result) {
						if(!list.contains(x.getId())) {
							list.add(x.getId());
							identifiers = identifiers + "," + x.getId();
						}
					}
					for(Xref x : result2) {
						if(!list.contains(x.getId())) {
							list.add(x.getId());
							if(x.getId().startsWith("MIMAT")) {
								mimat = x.getId();
							}
							identifiers = identifiers + "," + x.getId();
						}
					}
				} catch (IDMapperException e) {
					// could not be mapped
				}
			}
			identifiers = identifiers + "]";
			
			Node node = graph.addNode(miRNA);
			node.appendAttribute("identifiers", identifiers);
			node.appendAttribute("label", miRNA);
			node.appendAttribute("miRBaseAccession", mimat);
			node.appendAttribute("biologicalType", type);
			node.appendAttribute("organism", organism);
			countMiRNAs++;
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
