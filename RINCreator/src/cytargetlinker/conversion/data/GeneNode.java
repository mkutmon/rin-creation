package cytargetlinker.conversion.data;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.bridgedb.DataSource;
import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.Xref;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Node;

public class GeneNode {

	private String id;
	private List<String> entrez;
	private String label;
	private String biologicalType = "gene";
	
	private GeneNode() {
		entrez = new ArrayList<String>();
	}
	
	public static GeneNode createGeneNode(String id, String name, IDMapper mapper, DataSource in) {
		Xref ensembl;
		if(in.equals(DataSource.getBySystemCode("En"))) {
			ensembl = new Xref(id, in);
		} else {
			Xref xrefIn = new Xref(id, in);
			try {
				Set<Xref> ens = mapper.mapID(xrefIn, DataSource.getBySystemCode("En"));
				if(!ens.isEmpty()) {
					ensembl = ens.iterator().next();
				} else {
					return null;
				}
			} catch (IDMapperException e) {
				return null;
			}
		}
		
		GeneNode node = new GeneNode();
		node.setId(ensembl.getId());
		node.setLabel(name);
		try {
			Set<Xref> result = mapper.mapID(ensembl, DataSource.getBySystemCode("L"));
			for(Xref x : result) {
				String entrez = x.getId();
				if(!node.getEntrez().contains(entrez)) node.getEntrez().add(entrez);
			}
		} catch (IDMapperException e) {
		}
		
		return node;
	}

	public Node getNode(Graph graph) {
		Node node = graph.addNode(getId());
		node.appendAttribute("identifiers", getIdentifiers());
		node.appendAttribute("label", getLabel());
		node.appendAttribute("name", getLabel());
		node.appendAttribute("biologicalType", getBiologicalType());
		if(entrez.size() > 0) node.appendAttribute("entrezGeneID", entrez.get(0));
		return node;
	}
	
	private String getIdentifiers() {
		String ids = "[";
		ids = ids + getId();
		
		for(String name : getEntrez()) {
			ids = ids + "," + name;
		}
		
		ids = ids + "]";
		return ids;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public List<String> getEntrez() {
		return entrez;
	}

	public void setEntrez(List<String> entrez) {
		this.entrez = entrez;
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public String getBiologicalType() {
		return biologicalType;
	}

	public void setBiologicalType(String biologicalType) {
		this.biologicalType = biologicalType;
	}
}
