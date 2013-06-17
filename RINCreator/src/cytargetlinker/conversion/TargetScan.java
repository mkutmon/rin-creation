package cytargetlinker.conversion;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.bridgedb.DataSource;
import org.bridgedb.IDMapperStack;

import cytargetlinker.conversion.data.GeneNode;
import cytargetlinker.conversion.data.MTI;
import cytargetlinker.conversion.data.MiRNANode;
import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.utils.ArgsParser;
import cytargetlinker.conversion.utils.CommonAttributes;
import cytargetlinker.conversion.utils.Utils;
import cytargetlinker.conversion.utils.ArgsParser.AFilesAttributes;
import cytargetlinker.conversion.utils.ArgsParser.AFilesIn;
import cytargetlinker.conversion.utils.ArgsParser.AFilesOut;
import cytargetlinker.conversion.utils.ArgsParser.AHelp;
import cytargetlinker.conversion.utils.ArgsParser.GraphBuilder;

/**
 * Converts predicted miRNA target files from TargetScan to a XGMML network.
 * 
 * The input files for this script can be downloaded here:
 * human:
 * http://www.targetscan.org/vert_61/vert_61_data_download/Conserved_Site_Context_Scores.txt.zip
 * mouse:
 * http://www.targetscan.org/mmu_61/mmu_61_data_download/Conserved_Site_Context_Scores.txt.zip
 * 
 * @author martina
 */
public class TargetScan {

	private final static Logger log = Logger.getLogger(TargetScan.class.getName());
	private interface Args extends AHelp, AFilesIn, AFilesOut, AFilesAttributes {}
	private static Args pargs;
	
	/**
	 * MAIN METHOD
	 * USAGE: java -cp conversion.jar cytargetlinker.conversion.TargetScan 
	 * ARGUMENTS: 
	 * -i = input file
	 * -o = output file
	 * --bridgeDbFile = BridgeDb mapping files for source and target nodes
	 * --organism "Homo sapiens" or "Mus musculus" 
	 */
	public static void main(String[] args) throws Exception {
		
		pargs = ArgsParser.parse(args, Args.class);
		TargetScan converter = new TargetScan();
		converter.startConversion();

	}
	
	private IDMapperStack idMapper;
	private Graph graph;
	
	private Map<String, Integer> index;
	
	private Map<String, MiRNANode> miRNAs;
	private Map<String, GeneNode> genes;
	private Map<String, MTI> interactions;
	
	private List<String> nodesNotFound;
	private List<String> interactionIgnored;
	
	public TargetScan() throws Exception {
		index = new HashMap<String, Integer>();
		miRNAs = new HashMap<String, MiRNANode>();
		genes = new HashMap<String, GeneNode>();
		interactions = new HashMap<String, MTI>();
		nodesNotFound = new ArrayList<String>();
		interactionIgnored = new ArrayList<String>();
	}
	
	public void startConversion() throws Exception {
		File in = pargs.getInput();
		if(in != null) {
			Utils.setUpLogger(log, getLogFile(), false);

			if(pargs.isBridgeDbFiles()) idMapper = Utils.initIDMapper(pargs.getBridgeDbFiles(), false);
			if(idMapper == null) {
				log.severe("no bridgedb file specified");
			} else {
				if (pargs.getOrganism() != null) {
					log.info("conversion of TargetScan " + pargs.getOrganism() + " file started ...\n");
						
					ArgsParser.convertAndWrite(pargs, pargs, new GraphBuilder() {
						public Graph buildGraph(File in) throws Exception {
							return importTargetScan(in);
						}
					});
						
					log.info("conversion of TargetScan " + pargs.getOrganism() + " file finalized.\n");
				} else {
					log.severe("Please specify organism. --organism human or mouse.\n");
				}
			}
		} else {
			System.out.println("Please specify input file!");
		}
	}
	
	public Graph importTargetScan(File file) throws IOException {
		graph = new Graph();
		setNetworkAttributes(file);
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		String[] header = br.readLine().split("\t");
		
		index = new HashMap<String, Integer>();
		for(int i = 0; i < header.length; i++) index.put(header[i], i);
		
		log.info("Reading file start");
		String line = null;
		List<String[]> rows = new ArrayList<String[]>();
		while((line = br.readLine()) != null) {
			String [] str = line.split("\t", header.length);
			if(pargs.getOrganism().equals("Homo sapiens") && str[index.get("Gene Tax ID")].equals("9606")) {
				rows.add(line.split("\t", header.length));
			} else if(pargs.getOrganism().equals("Mus musculus") && str[index.get("Gene Tax ID")].equals("10090")) {
				rows.add(line.split("\t", header.length));
			}
		}
		log.info("Reading file finsihed");
		
		log.info("Create interactions start");
		int countNotMapped = 0;
		for(String[] row : rows) {
			String geneId = createGeneNode(row);
			String miRNA = createMiRNANode(row);
			if(geneId != null && miRNA != null) {
				if(!interactions.containsKey(miRNA + "_" + geneId)) {
					String score = row[index.get("context+ score")];
					try{
						Double.valueOf(score);
					} catch(NumberFormatException e) {
						score = "";
					}
					MTI mti = new MTI(miRNAs.get(miRNA), genes.get(geneId), score);
					interactions.put(miRNA + "_" + geneId, mti);
				}
			} else {
				if(!interactionIgnored.contains(row[index.get("miRNA")] + " -> " + row[index.get("Gene ID")])) {
					interactionIgnored.add(row[index.get("miRNA")] + " -> " + row[index.get("Gene ID")]);
				}
				countNotMapped++;
			}
		}
		log.info("Create interactions finished (" + countNotMapped + " were not created)");
		
		int count = 0;
		for(String str : interactions.keySet()) {
			MTI mti = interactions.get(str);
			mti.createEdge(graph, "TargetScan version 6.2", count);
			count++;
		}
		
		log.info(interactions.size() + " interactions have been found.\n" + genes.size() + " gene nodes.\n" + miRNAs.size() + " miRNA nodes.\n");
		log.info(interactionIgnored.size() + " interactions were ignores because " + nodesNotFound.size() + " nodes could not be mapped");
		
		cleanUp();
		return graph;
	}
	
	private String createGeneNode(String[] r) {
		String geneName = r[index.get("Gene Symbol")];
		String geneId = r[index.get("Gene ID")];

		GeneNode node = GeneNode.createGeneNode(geneId, geneName, idMapper, DataSource.getBySystemCode("L"));
		if(node != null) {
			if(!genes.containsKey(node.getId())) {
				genes.put(node.getId(), node);
			}
			return node.getId();
		} else {
			if(!nodesNotFound.contains(r[index.get("Gene ID")])) {
				nodesNotFound.add(r[index.get("Gene ID")]);
			}
			return null;
		}
	}
	
	
	
	private String createMiRNANode(String[] r) {
		
		MiRNANode node = MiRNANode.createMiRNANode(r[index.get("miRNA")], idMapper, DataSource.getBySystemCode("Mb"));
		if(node != null) {
			if(!miRNAs.containsKey(node.getId())) {
				miRNAs.put(node.getId(), node);
			}
			return node.getId();
		} else {
			if(!nodesNotFound.contains(r[index.get("miRNA")])) {
				nodesNotFound.add(r[index.get("miRNA")]);
			}
			return null;
		}
	}
	
	
	private void cleanUp() {
		genes.clear();
		miRNAs.clear();
		interactions.clear();
		nodesNotFound.clear();
		interactionIgnored.clear();
	}
	
	private void setNetworkAttributes(File input) {
		graph.setTitle("TargetScan " + pargs.getOrganism() + " version 6.2");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), "TargetScan " + pargs.getOrganism() + " version 6.2");
		graph.setAttribute(CommonAttributes.TYPE.getName(), "predicted miRNA targets");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), input.getName());
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), "EntrezGene");
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "miRBase");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "mirna");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "gene");
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
