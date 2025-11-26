import gdown
import os
import re
import requests

def get_file_ids_from_folder(folder_id):
    """
    Scrape file IDs from a public Google Drive folder.
    """
    files = []
    
    # Access the folder page
    url = f"https://drive.google.com/drive/folders/{folder_id}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        content = response.text
        
        # Extract file IDs and names using regex
        # Pattern matches the data structure Google Drive uses
        pattern = r'\["([^"]{33}|[^"]{19})","([^"]+)","application/[^"]*"\]|\["([^"]{33}|[^"]{19})","([^"]+)","[^a][^p][^p][^"]*"\]'
        matches = re.findall(pattern, content)
        
        for match in matches:
            file_id = match[0] or match[2]
            file_name = match[1] or match[3]
            if file_id and file_name and len(file_id) >= 19:
                files.append({'id': file_id, 'name': file_name})
        
        # Remove duplicates
        seen = set()
        unique_files = []
        for f in files:
            if f['id'] not in seen:
                seen.add(f['id'])
                unique_files.append(f)
        
        return unique_files
        
    except Exception as e:
        print(f"Error scraping folder: {e}")
        return []

def download_folder_files(folder_url, output_dir="downloads", skip_existing=True):
    """
    Download all files from a Google Drive folder.
    
    Args:
        folder_url: The Google Drive folder URL or folder ID
        output_dir: Directory where files will be downloaded
        skip_existing: Skip files that already exist locally
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract folder ID from URL if needed
    if "folders/" in folder_url:
        folder_id = folder_url.split("folders/")[1].split("?")[0]
    else:
        folder_id = folder_url
    
    print(f"Fetching file list from folder: {folder_id}")
    print(f"Output directory: {output_dir}\n")
    
    # Get list of files
    files = get_file_ids_from_folder(folder_id)
    
    if not files:
        print("⚠️  No files found using scraping method.")
        print("Trying alternative method with gdown...\n")
        
        # Fallback: try to use gdown's download with the folder directly
        # Download directly to the specified output directory
        try:
            gdown.download_folder(id=folder_id, output=output_dir, quiet=False, remaining_ok=True)
            print(f"\n✅ Files downloaded to: {output_dir}")
            return
        except Exception as e:
            print(f"❌ Alternative method also failed: {e}")
            return
    
    print(f"Found {len(files)} files\n")
    
    # Show file list
    print("Files to download:")
    for i, f in enumerate(files, 1):
        print(f"  {i}. {f['name']}")
    print()
    
    # Download each file
    success_count = 0
    skip_count = 0
    fail_count = 0
    
    for i, file_info in enumerate(files, 1):
        file_id = file_info['id']
        file_name = file_info['name']
        output_path = os.path.join(output_dir, file_name)
        
        # Skip if file already exists
        if skip_existing and os.path.exists(output_path):
            print(f"[{i}/{len(files)}] ⏭️  Skipping (exists): {file_name}")
            skip_count += 1
            continue
        
        print(f"[{i}/{len(files)}] Downloading: {file_name}")
        
        try:
            gdown.download(id=file_id, output=output_path, quiet=False)
            success_count += 1
        except Exception as e:
            print(f"  ❌ Error: {e}")
            fail_count += 1
            continue
        
        print()
    
    print(f"\n{'='*60}")
    print(f"✅ Download complete!")
    print(f"   Success: {success_count}")
    print(f"   Skipped: {skip_count}")
    print(f"   Failed:  {fail_count}")
    print(f"{'='*60}")

if __name__ == "__main__":
    # Replace with your Google Drive folder URL or ID
    FOLDER_URL = "https://drive.google.com/drive/folders/1014qkKPpxFi5oHS1GjPOTf8gm8fJN839"
    
    # Specify output directory
    OUTPUT_DIR = "/biodata/dep_mercier/grp_marques/marques/Carnivorous-plants/RNA_seq/D_capensis"
    
    download_folder_files(FOLDER_URL, OUTPUT_DIR)
