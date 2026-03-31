import pandas as pd
import glob
import os
import argparse

def build_html():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max_price", type=float, default=999999.0)
    args = parser.parse_args()

    # 1. Extract Scores
    print("Extracting docking scores from results folder...")
    zinc_scores = {}
    pdbqt_files = glob.glob(os.path.join("results", "**", "*.pdbqt"), recursive=True)
    
    for file in pdbqt_files:
        best_score = None
        zinc_id = None
        try:
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith("REMARK VINA RESULT:") and best_score is None:
                        best_score = float(line.split()[3])
                    elif line.startswith("REMARK  Name ="):
                        zinc_id = line.split("=")[1].strip()
                        break 
            if zinc_id and best_score is not None:
                zinc_scores[zinc_id] = best_score
        except Exception:
            pass

    # 2. Load Pricing
    try:
        df = pd.read_csv("MASTER_PRICING_SORTED.csv")
    except FileNotFoundError:
        print("Run merge_pricing.py first!")
        return
        
    print(f"Building dashboard...")
    
    html_header = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Ligand Dashboard</title>
        <script src="https://unpkg.com/smiles-drawer@1.0.10/dist/smiles-drawer.min.js"></script>
        <style>
            body { font-family: 'Segoe UI', sans-serif; background: #f4f7f6; padding: 20px; }
            .grid-container { display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 20px; }
            .card { background: white; border-radius: 8px; padding: 15px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); 
                    display: flex; flex-direction: column; align-items: center; }
            canvas { width: 250px; height: 200px; }
            .info { width: 100%; text-align: left; margin-top: 10px; border-top: 1px solid #eee; padding-top: 10px;}
            .title { font-size: 1.2em; font-weight: bold; color: #1a237e; margin-bottom: 5px; text-align: center; }
            .score { font-size: 1.1em; font-weight: bold; color: #d32f2f; }
            .price { font-size: 1.1em; font-weight: bold; color: #2e7d32; }
            .meta { font-size: 0.85em; color: #666; margin-top: 5px; }
        </style>
    </head>
    <body>
        <div class="grid-container">
    """

    cards_html = ""

    for idx, row in df.iterrows():
        ppg = row.get('Price_Per_Gram')
        if pd.isna(ppg) or ppg > args.max_price: continue
            
        smiles = row.get('Ligand_SMILES', '')
        orig_id = row.get('Original_ID', 'Unknown')
        vend_id = row.get('Vendor_ID', 'Unknown')
        
        # Determine the display title
        display_title = orig_id if orig_id != 'Unknown ID' else vend_id
        
        # Link the score using the Original ZINC ID
        score = zinc_scores.get(orig_id)
        score_text = f"{score:.1f}" if score is not None else "N/A"
        
        card = f"""
            <div class="card">
                <div class="title">{display_title}</div>
                <canvas data-smiles="{smiles}"></canvas>
                <div class="info">
                    <div class="score">Score: {score_text}</div>
                    <div class="price">Price: ${ppg:.2f} / g</div>
                    <div class="meta">
                        <b>Vendor:</b> {row.get('Vendor_Name', 'Unknown')}<br>
                        <b>Pack:</b> ${row.get('Price_USD', 0):.2f} for {row.get('Pack_Size', 'Unknown')}<br>
                        <b>Catalog ID:</b> {vend_id}
                    </div>
                </div>
            </div>
        """
        cards_html += card

    html_footer = """
        </div>
        <script>
            let options = { width: 250, height: 200, terminalCarbons: true };
            let smilesDrawer = new SmilesDrawer.Drawer(options);
            let canvases = document.querySelectorAll('canvas[data-smiles]');
            canvases.forEach(function(canvas) {
                let sm = canvas.getAttribute('data-smiles');
                if (sm) {
                    SmilesDrawer.parse(sm, function(tree) {
                        smilesDrawer.draw(tree, canvas, 'light', false);
                    });
                }
            });
        </script>
    </body>
    </html>
    """

    with open("dashboard.html", "w") as f:
        f.write(html_header + cards_html + html_footer)

    print("Success! Dashboard generated.")

if __name__ == "__main__":
    build_html()