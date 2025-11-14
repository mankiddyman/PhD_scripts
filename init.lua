-- Bootstrap lazy.nvim
local lazypath = vim.fn.stdpath("data") .. "/lazy/lazy.nvim"
if not (vim.uv or vim.loop).fs_stat(lazypath) then
  vim.fn.system({
    "git",
    "clone",
    "--filter=blob:none",
    "https://github.com/folke/lazy.nvim.git",
    "--branch=stable",
    lazypath,
  })
end
vim.opt.rtp:prepend(lazypath)

-- Plugins
local opts = {}
local plugins = {
  { "catppuccin/nvim", name = "catppuccin", priority = 1000 },
  { "nvim-treesitter/nvim-treesitter", build = ":TSUpdate" },
  { "github/copilot.vim", config = function() end },
  { "catgoose/nvim-colorizer.lua", event = "BufReadPre", opts = {} },

  -- vim-slime
	 {'jpalardy/vim-slime',
          config=function()
		  vim.cmd([[let g:slime_target="tmux"]])
			vim.cmd([[let g:slime_paste_file="$HOME/.slime_paste"]])
			vim.cmd([[let g:slime_target_pane="1"]])
	  end,
	  },
}

require("lazy").setup(plugins, opts)

-- Basic settings
vim.g.mapleader = " "
vim.opt.expandtab = true
vim.opt.tabstop = 2
vim.opt.softtabstop = 2
vim.opt.shiftwidth = 2
vim.opt.number = true
vim.opt.relativenumber = true

-- Catppuccin setup
require("catppuccin").setup()
vim.cmd.colorscheme("catppuccin")

-- Treesitter setup
require("nvim-treesitter.configs").setup({
  ensure_installed = { "lua", "bash", "python", "r" },
  highlight = { enable = true },
  indent = { enable = true },
})

-- Keymaps
-- Save
vim.keymap.set('n', '<leader>w', ':w<CR>', { desc = 'Save File' })
vim.keymap.set('v', '<leader>w', ':w<CR>', { desc = 'Save File' })

-- Select all
vim.keymap.set('n', '<leader>a', 'ggVG', { desc = 'Select All' })
vim.keymap.set('v', '<leader>a', 'ggVG', { desc = 'Select All' })

-- Go to end/start of line
vim.keymap.set('n', '<leader>e', '$', { desc = 'Go to End of Line' })
vim.keymap.set('n', '<leader>b', '0', { desc = 'Go to Beginning of Line' })

-- OSC 52 clipboard support for SSH
local function copy_to_clipboard(lines)
  local content = table.concat(lines, '\n')
  local base64_content = vim.base64.encode(content)
  local osc52 = '\x1b]52;c;' .. base64_content .. '\x1b\\'
  io.stdout:write(osc52)
  io.stdout:flush()
  vim.notify('Copied to clipboard via OSC 52', vim.log.levels.INFO)
  return content
end

vim.g.clipboard = {
  name = 'OSC 52',
  copy = {
    ['+'] = copy_to_clipboard,
    ['*'] = copy_to_clipboard,
  },
  paste = {
    ['+'] = function() return {vim.fn.getreg('"'), vim.fn.getregtype('"')} end,
    ['*'] = function() return {vim.fn.getreg('"'), vim.fn.getregtype('"')} end,
  },
}

-- Map Ctrl+C to yank (copy)
vim.keymap.set('v', '<C-c>', '"+y', { desc = 'Copy to Clipboard' })
vim.keymap.set('n', '<C-c>', '"+yy', { desc = 'Copy line to Clipboard' })
-- also for visual line mode

-- Map Ctrl+V to paste
vim.keymap.set('n', '<C-v>', 'p', { desc = 'Paste' })
vim.keymap.set('i', '<C-v>', '<C-r>"', { desc = 'Paste' })



-- insert current date in format like this 12:50|Fr.|Okt.|31|2025 using leader D
vim.keymap.set('n', '<leader>d', function()
  local date_str = os.date("%H:%M|%a.|%b.|%d|%Y")
  vim.api.nvim_put({date_str}, 'c', true, true)
end, { desc = 'Insert Current Date' })

