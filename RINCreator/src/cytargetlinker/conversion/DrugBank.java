package cytargetlinker.conversion;

import java.io.File;
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
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.Namespace;
import org.jdom.input.SAXBuilder;

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
 * Converts drug target xml file from DrugBank to a XGMML network.
 * 
 * The input files for this script can be downloaded here:
 * http://www.drugbank.ca/downloads
 * 
 * @author Tina
 */
public class DrugBank {
	private final static Logger log = Logger.getLogger(MiRecords.class.getName());
	private static Args pargs;
	private interface Args extends AHelp, AFilesIn, AFilesOut, AFilesAttributes {}
	
	/**
	 * USAGE: java -cp conversion.jar cytargetlinker.conversion.DrugBank 
	 * ARGUMENTS: 
	 * -i = input file
	 * [-o] = output file
	 * [--bridgeDbFile] = BridgeDb mapping files for source and target nodes
	 * [-l] = log file
	 * [-t] = type (xgmml or gml)
	 */
	public static void main(String argv[]) throws Exception {
		pargs = ArgsParser.parse(argv, Args.class);
		DrugBank converter = new DrugBank();
		converter.startConversion();
	}

	private IDMapperStack idMapper;
	private Graph graph;
	private Map<String, List<String>> edges;
	private List<String> foundConnections;
	
	private int countEdge = 0;
	private int countTargets = 0;
	private int countDrugs = 0;
	private boolean mapping = false;
	
	private Map<String, Drug> drugs;
	private Map<String, Target> targets;
	private Namespace nsDrugBank = Namespace.getNamespace("http://drugbank.ca");

	public DrugBank() throws Exception {
		edges = new HashMap<String, List<String>>();
		foundConnections = new ArrayList<String>();
		drugs = new HashMap<String, DrugBank.Drug>();
		targets = new HashMap<String, DrugBank.Target>();
	}
	
	private void startConversion() throws Exception {
		File in = pargs.getInput();
		if(in != null && in.getName().endsWith(".xml")) {
			Utils.setUpLogger(log, getLogFile(), false);
			
			if(pargs.isBridgeDbFiles()) {
				idMapper = Utils.initIDMapper(pargs.getBridgeDbFiles(), false);
				if(idMapper != null) {
					mapping = true;
					log.info("identifier mapping on");
				}
			}
			
			log.info("conversion of drugbank started ...\n");
			
			ArgsParser.convertAndWrite(pargs, pargs, new GraphBuilder() {
				public Graph buildGraph(File in) throws Exception {
					return convert(in);
				}
			});
				
			log.info("Drugbank network file created.\n\n");
			
		} else {
			log.severe("Please specify XML input file!");
		}
	}
	
	
	
	@SuppressWarnings("unchecked")
	private Graph convert(File in) throws JDOMException, IOException {
		graph = new Graph();
		
		SAXBuilder builder = new SAXBuilder();
		try { 
			Document document = (Document) builder.build(in);
			Element rootNode = document.getRootElement();
			
			List<Element> list = rootNode.getChildren("drug", nsDrugBank);
			for(Element e : list) {
				
				String dbId = e.getChild("drugbank-id", nsDrugBank).getValue();
				String name = e.getChild("name", nsDrugBank).getValue();
				String cas = e.getChild("cas-number", nsDrugBank).getValue();

				Drug drug = new Drug(dbId, name, cas);
				getXrefs(e, drug);
				Element targets = e.getChild("targets", nsDrugBank);
				List<Element> targetList = targets.getChildren("target", nsDrugBank);
				
				for(Element t : targetList) {
					String partner = t.getAttributeValue("partner");
					
					List<String> listActions = new ArrayList<String>();
					Element actions = t.getChild("actions", nsDrugBank);
					List<Element> actionList = actions.getChildren("action", nsDrugBank);
					for(Element a : actionList) {
						String action = a.getValue();
						listActions.add(action);
					}
					
					drug.addTargetId(partner, listActions);
				}
				
				drugs.put(dbId, drug);
			}

			Element partners = rootNode.getChild("partners", nsDrugBank);
			List<Element> partnerList = partners.getChildren("partner", nsDrugBank);

			for(Element e : partnerList) {
				String geneName = e.getChild("gene-name", nsDrugBank).getValue();
				String id = e.getAttributeValue("id");
				
				Target target = new Target(geneName, id);
				
				// if a target has a uniprot id
				if(getXrefs(e, target)) {
					targets.put(id, target);
				}
			}
		} catch (IOException io) {
			 System.out.println(io.getMessage());
		 }
		
		for(String dbId : drugs.keySet()) {
			Drug drug = drugs.get(dbId);
			for(String id : drug.getTargetIds().keySet()) {
				Target target = targets.get(id);
				// might be null if the target does not have any xrefs
				if(target != null) {
					drug.getTargets().put(target, drug.getTargetIds().get(id));
				}
			}
		}
		
		setNetworkAttributes(in);
		
		for(String dbId : drugs.keySet()) {
			createDrugNode(drugs.get(dbId));
			for(Target target : drugs.get(dbId).getTargets().keySet()) {
				createTargetNode(target);
				createEdge(dbId, target.getId());
			}
		}
		
		log.info(foundConnections.size() + " interactions have been found.\n" + countDrugs + " drugs.\n" + countTargets + " target nodes.\n");
		return graph;
	}

	private void createDrugNode(Drug drug) {
		if(graph.getNode(drug.getDbId()) == null) {
			String identifiers = "[" + drug.getDbId();
			
			if(drug.getCas() != null && !drug.getCas().equals("")) {
				identifiers = identifiers + "," + drug.getCas();
			}

			identifiers = identifiers + "]";
			
			Node node = graph.addNode(drug.getDbId());
			node.appendAttribute("identifiers", identifiers);
			node.appendAttribute("label", drug.getDbId());
			node.appendAttribute("name", drug.getName());
			node.appendAttribute("biologicalType", "drug");
			node.appendAttribute("cas-number", drug.getCas());
			
			for(String str : drug.getRefs().keySet()) {
				node.appendAttribute(str, drug.getRefs().get(str));
			}
			
			countDrugs++;
		}
	}
	
	private String createTargetNode(Target target) {
		String id = target.getRefs().get("UniProtKB");
		
		if(id != null && graph.getNode(id) == null) {
			String identifiers = "[";
			identifiers = identifiers + id;

			String ensembl = "";
			String entrez = "";
			if(mapping) {
				Xref xrefIn = new Xref(id, DataSource.getBySystemCode("S"));
				try {
					Set<Xref> e = idMapper.mapID(xrefIn, DataSource.getBySystemCode("En"));
					if(!e.isEmpty()) {
						ensembl = e.iterator().next().getId();
						identifiers = identifiers + "," + ensembl;
					}
					Set<Xref> result = idMapper.mapID(xrefIn, DataSource.getBySystemCode("L"));
					if(!result.isEmpty()) {
						entrez = result.iterator().next().getId();
						identifiers = identifiers + "," + entrez;
					}
				} catch (IDMapperException ex) {
					// could not be mapped
				}
			}
			
			identifiers = identifiers + "]";
			
			Node node = graph.addNode(id);
			node.appendAttribute("identifiers", identifiers);
			node.appendAttribute("label", target.getName());
			node.appendAttribute("name", target.getName());
			node.appendAttribute("ensemblID", ensembl);
			node.appendAttribute("entrezGeneID", entrez);
			node.appendAttribute("uniprot", id);
			node.appendAttribute("biologicalType", "drug target");
		
			for(String str : target.getRefs().keySet()) {
				node.appendAttribute(str, target.getRefs().get(str));
			}
			
			countTargets++;
		}
		return id;
	}
	
	private void setNetworkAttributes(File in) {
		graph.setTitle("DrugBank");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), "DrugBank");
		graph.setAttribute(CommonAttributes.TYPE.getName(), "drug-target interactions");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), in.getName());
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), "UniProt/TrEMBL");
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "DrugBankID");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "drug");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "drug target");
	}

	private void createEdge(String dbId, String id) {
		if(edges.containsKey(dbId)) {
			if(!edges.get(dbId).contains(id)) {
				Edge e = graph.addEdge("" + countEdge, graph.getNode(dbId),graph.getNode(id));
				e.setAttribute("datasource", "DrugBank");
				e.setAttribute("interactionType", "drug-target");
				
				edges.get(dbId).add(id);
				foundConnections.add(dbId + "\t" + id);
				countEdge++;
			}
		} else {
			Edge e = graph.addEdge("" + countEdge, graph.getNode(dbId),graph.getNode(id));
			e.setAttribute("datasource", "DrugBank");
			e.setAttribute("interactionType", "drug-target");
			
			List<String> list = new ArrayList<String>();
			list.add(id);
			foundConnections.add(dbId + "\t" + id);
			edges.put(dbId, list);
			countEdge++;
		}
	}

	@SuppressWarnings("unchecked")
	private void getXrefs(Element drugElement, Drug drug) {
		Element refs = drugElement.getChild("external-identifiers", Namespace.getNamespace("http://drugbank.ca"));
		List<Element> refList = refs.getChildren("external-identifier", Namespace.getNamespace("http://drugbank.ca"));
		for(Element e : refList) {
			drug.getRefs().put(e.getChild("resource", Namespace.getNamespace("http://drugbank.ca")).getValue(), e.getChild("identifier", Namespace.getNamespace("http://drugbank.ca")).getValue());
		}
	}
	
	@SuppressWarnings("unchecked")
	private boolean getXrefs(Element drugElement, Target target) {
		Element refs = drugElement.getChild("external-identifiers", Namespace.getNamespace("http://drugbank.ca"));
		List<Element> refList = refs.getChildren("external-identifier", Namespace.getNamespace("http://drugbank.ca"));
		for(Element e : refList) {
			target.getRefs().put(e.getChild("resource", Namespace.getNamespace("http://drugbank.ca")).getValue(), e.getChild("identifier", Namespace.getNamespace("http://drugbank.ca")).getValue());
		}
		if(target.getRefs().containsKey("UniProtKB")) {
			target.setId(target.getRefs().get("UniProtKB"));
			return true;
		}
		return false;
	}

	/**
	 * class that stores information about a target
	 */
	public class Target {
		private String name;
		private String id;
		private String partnerId;
		private Map<String, String> refs;
		
		public Target(String name, String id) {
			super();
			this.name = name;
			this.partnerId = id;
			refs = new HashMap<String, String>();
		}

		public void setId(String id) {
			this.id = id;
		}

		public String getName() {
			return name;
		}

		public String getId() {
			return id;
		}

		public Map<String, String> getRefs() {
			return refs;
		}

		public String getPartnerId() {
			return partnerId;
		}
	}
	
	/**
	 * class that stores information about a drug
	 */
	private class Drug {
		private String dbId;
		private String name;
		private String cas;
		private Map<String, String> refs;
		private Map<String, List<String>> targetIds;
		private Map<Target, List<String>> targets;

		public Drug(String dbId, String name, String cas) {
			super();
			this.dbId = dbId;
			this.name = name;
			this.cas = cas;
			targetIds = new HashMap<String, List<String>>();
			refs = new HashMap<String, String>();
			targets = new HashMap<DrugBank.Target, List<String>>();
		}
		
		public void addTargetId(String id, List<String> actions) {
			if(!targetIds.containsKey(id)) {
				targetIds.put(id, actions);
			}
		}

		public String getDbId() {
			return dbId;
		}

		public String getName() {
			name = name.replace("[", "");
			name = name.replace("]", "");
			return name;
		}

		public String getCas() {
			return cas;
		}

		public Map<String, List<String>> getTargetIds() {
			return targetIds;
		}

		public Map<String, String> getRefs() {
			return refs;
		}

		public Map<Target, List<String>> getTargets() {
			return targets;
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