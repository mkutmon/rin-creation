package cytargetlinker.conversion.utils;
import java.util.HashMap;
import java.util.Map;

/**
 * adds species code if it is missing in the miRNA identifier
 * @author Thomas
 *
 */
public class SpeciesCodes {
	private static Map<String, String> species2code = new HashMap<String, String>();
	static {
		//TarBase
		species2code.put("C. elegans", "cel");
		species2code.put("D. rerio", "dre");
		species2code.put("Drosophila", "dme");
		species2code.put("Human", "hsa");
		species2code.put("Mouse", "mmu");
		species2code.put("Rat", "rno");
		
		//MiRecords
		species2code.put("Homo sapiens", "hsa");
		species2code.put("Caenorhabditis elegans", "cel");
		species2code.put("Drosophila melanogaster", "dme");
		species2code.put("Mus musculus", "mmu");
		species2code.put("Danio rerio", "dre");
		species2code.put("Ovis aries", "oar");
		species2code.put("Rattus norvegicus", "rno");
		species2code.put("Mus musculus miR-24", "mmu");
		species2code.put("Kaposi sarcoma-associated herpesvirus", "kshv");
		species2code.put("Bombyx mori", "bmo");
		species2code.put("Bos taurus", "bta");
		species2code.put("Xenopus laevis", "xla");
		species2code.put("Gallus gallus", "gga");
		species2code.put("Trichomonas vaginalis", "tva");
		
		// miRDB
		species2code.put("Canis lupus familiaris", "cfa");
	}
	
	public static String getCode(String species) {
		return species2code.get(species);
	}
}
